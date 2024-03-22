
# This script can be used to replicate the results of the paper 'How to measure
# publication relatedness using MeSH-terms', which will be submitted for review
# to the 28th International Conference on Science, Technology and Innovation
# Indicators, 2024.

# In most parts of the code, the global environment is not cleaned, i.e.,
# intermediate objects are not removed. If RAM is needed, the global environment
# can be cleaned before running each chapter (indicated by 5 consecutive lines
# starting with '#), since all objects necessary to execute a chapter are loaded
# within that chapter.

# Author: Emil Dolmer Alnor, Aarhus University, ea@ps.au.dk

library(httr)
library(jsonlite)
library(rentrez)
library(xml2)
library(XML)
library(stringr)
library(igraph)
library(tidyr)
library(dplyr)
library(matrixStats)
library(Matrix)
library(future.apply)
library(proxyC)
library(effsize)
library(patchwork)
library(ggplot2)
library(viridis)
library(openxlsx)
library(readxl)

# ***************************************************************** #
#*******************************************************************#
### MeSH Hiearchy                                                   #
#*******************************************************************#
# ***************************************************************** #

#Download data
download.file(
  "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/d2024.bin",
  destfile = "d2024.bin",
  mode = "wb"
)

lines <- readLines("d2024.bin")

#Mesh + Muid
mh <- lines[str_detect(lines, "MH = ")] %>% str_remove("MH = ")
muid <- lines[str_detect(lines, "UI = ")] %>% str_remove("UI = ")

#Tree numbers
nrl <- which(lines == "*NEWRECORD")
tnl <- list()

for (i in 1:length(nrl)) {
  j <- nrl[i]
  
  if (j == max(nrl)) k <- nrl[length(nrl)]
  else k <- nrl[i+1]
  
  s <- lines[j:k]
  
  t <- list(s[str_detect(s, "MN = ")]) %>%
    lapply(function(x) gsub("MN = ", "", x))
  
  tnl[[i]] <- t
}

tn <- sapply(
  tnl,
  function(x) unlist(x) %>% paste(collapse = ";")
) 

tree <- data.frame(muid, mh, tn) %>% filter(tn != "")

#Descendant
tree$pattern <- paste0(
  str_replace_all(tree$tn, ";", "\\\\.|"),
  "\\."
)

fn <- function(x) {
  
  index <- which(str_detect(tree$tn, x))
  
  nodes <- paste(
    tree$muid[index],
    collapse = ";"
  ) %>%
    ifelse(. == "", NA, .)
}

plan(multisession, workers = 6) 
tree$desc <- future_sapply(tree$pattern, fn)

#Children
tree <- tree %>%
  mutate(
    tnc1 = paste0(tn, "\\.\\d{1,5}$") %>%
      str_replace_all(";", "\\\\.\\\\d{1,5}$|"),
    tnc2 = str_replace_all(tn, ';', "\\\\.\\\\d{1,5};|") %>%
      paste0("\\.\\d{1,5};"),
    tnc = paste(tnc1, tnc2, sep = '|')
  )

tree$chld <- future_sapply(tree$tnc, fn)

tree <- tree %>% select(muid, mh, tn, desc, chld)

save(tree, file = "tree.rda")

# ********************************************* #
#  Edgelist                                     #
# ********************************************* #

load("tree.rda")

#We start by creating a data.frame with one row for each node representation of the MeSH-terms.
nodes <- tree %>% separate_rows(tn, sep = ';') %>% select(muid, tn)

#Next we modify this data.frame to show the tree number of the parent for each node. Apart from highest lvl MeSH-terms (just below the categories) the parent is their tree number with the last digits and '.' removed
children <- nodes %>% 
  mutate(
    tnp = ifelse(
      str_detect(tn, '\\.'), 
      str_remove(tn, '\\.\\d+$'), 
      str_extract(tn, '^.')
    )
  ) %>% select(muid, tnp) %>% rename(chld = muid)

#We can now add the children to each of the MeSh-terms.
edgelist <- nodes %>% 
  inner_join(children, by = c('tn' = 'tnp')) %>% 
  distinct(chld, muid, .keep_all = T) #Many parent-child-pairs are parent-child-pairs with more than 1 node.

#We now need to add that the parents of the highest level MeSH-terms are the categories 
children <- edgelist %>% 
  filter(!str_detect(tn, '\\.')) %>% #All nodes with '.' are not highest lvl
  distinct(muid, tn) %>%  
  mutate(tn = str_extract(tn, '^.')) %>%  #We take the category letter
  rename(chld = muid) %>% 
  mutate(muid = tn)

edgelist <- edgelist %>% rbind(children)

#Finally we add that the parent of the categories is the tree
root <- data.frame(
  chld = unique(children$tn),
  tn = 'ROOT',
  muid = 'ROOT'
)

edgelist <- edgelist %>% rbind(root)

save(edgelist, file = 'edgelist.rda')

# ***************************************************************** #
#*******************************************************************#
### TREC Genomics 06'                                               #
#*******************************************************************#
# ***************************************************************** #

#topics
download.file(
  "https://dmice.ohsu.edu/trec-gen/data/2006/topics/2006topics.xls",
  destfile = "2006topics.xls",
  mode = "wb"
)

temp <- read_excel("2006topics.xls")

temp_fn <- function(type, first, last) {
  
  mat <- temp %>%
    slice(first:last) %>%
    mutate(type = type)
  
  colnames(mat) <- c('nid', 'id', 'gene', 'facet2', 'need', 'type')
  
  return(mat)
}

topics6 <- rbind(
  temp_fn("disease", 3, 9),
  temp_fn("process", 14, 20),
  temp_fn("function", 25, 31),
  temp_fn("impact", 36, 42)
) %>%
  select('nid', 'id', 'type', 'need', 'gene', 'facet2')

write.xlsx(topics6, "topics6.xlsx")

t6 <- read_xlsx("topics6.xlsx")

#rjs
download.file(
  "https://dmice.ohsu.edu/trec-gen/data/2006/trec2006.raw.relevance.tsv.txt",
  destfile = "trec2006rrjs.txt",
  mode = "wb"
)

rjs6 <- read.delim("trec2006rrjs.txt", header = F) %>%
  slice(4:nrow(.)) %>%
  select(topic = V1, pmid = V2, rj = V6) %>%
  mutate(rj = recode(rj, "NOT" = 0, "POSSIBLY" = 1, "DEFINITELY" = 2)) %>% 
  group_by(pmid, topic) %>% 
  summarise(rj = max(rj))

write.csv(rjs6, "rjs6.txt", row.names = F)

rjStats <- rjs6 %>%
  group_by(topic) %>%
  summarise(
    no       = sum(rj ==0 ),
    possibly = sum(rj == 1),
    yes      = sum(rj == 2)
  ) %>%
  mutate(noRat = round(no / (no + possibly + yes), 3)) %>% 
  arrange(noRat)

save(rjStats, file = 'rjStats.rda')

#pmids
download.file(
  "https://dmice.ohsu.edu/trec-gen/data/2006/medline/pmids.txt",
  destfile = "pmids6.txt",
  mode = "wb"
)

pmids6 <- read.delim("pmids6.txt", header = T)

write.csv(pmids6, "pmids6.txt", row.names = F)

pmids6 <- read.csv("pmids6.txt") %>% unlist()

download.file(
  "https://dmice.ohsu.edu/trec-gen/data/2006/medline/error_pmids.txt",
  destfile = "error_pmids.txt",
  mode = "wb"
)

error_pmids <- read.delim("error_pmids.txt", header = F) %>%
  unlist()

pmids6_clean <- setdiff(pmids6, error_pmids)

write.csv(pmids6_clean, "pmids6_clean.txt", row.names = F)

# ***************************************************************** #
#*******************************************************************#
### Publication MeSH data                                           #
#*******************************************************************#
# ***************************************************************** #

# An API-key for Entrez is required in this section. The API-key is used to
# speed up the download. To replicate the results without registering for an
# API-key, remove the line where 'ak' is defined and remove', api_key = ak' from
# the 'entrez_fetch'.

# ********************************************* #
#  Functions                                    #
# ********************************************* #

# Sleep
sleep <- function(time = 1) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (elapsed < time) {
    Sys.sleep(time - elapsed)
  }
}

# Progress bar
progress <- function(iterator, sequence) {
  percentage <- (which(sequence == iterator)-1) / length(sequence) * 100
  timestamp <- format(Sys.time(), format = "%d-%b %H:%M")
  cat("Round:", which(sequence == iterator), "/", length(sequence),
      "   Progress:", round(percentage, 0), "%",
      "  ", timestamp, "\n")
}

# Extract
fnQ <- function(target) {
  result <- sapply(
    mh,
    \(x) {
      if (length(x) == 1) NA
      else if (length(x) == 2) x[['QualifierName']][['.attrs']][[target]]
      else {
        elements <- x[-1] #Remove descriptor from list
        UIs <- sapply(elements, \(y) y[['.attrs']][[target]]) %>% 
          paste(collapse = ';')
      }
    }
  ) 
}

fnM <- function(target) {
  result <- sapply(mh, \(x) x[['DescriptorName']][['.attrs']][[target]])
}

# ********************************************* #
#  Data and loop                                #
# ********************************************* #
results <- list()
min <- 1
max <- 50
j <- 1
pmids <- read.csv("pmids6_clean.txt") %>% unlist()
ak <- readLines("C:/Users/au544242/OneDrive - Aarhus universitet/Kodning/rentrez_token.txt")

while (T == T) {
  
  progress(min, 1:length(pmids))
  start_time <- Sys.time()
  
  pmid <- pmids[min:max]
  
  xml <- entrez_fetch(
    db = "pubmed", id = pmid, rettype = "xml", parsed = T, api_key = ak) 
  
  list <- xmlToList(xml) 
  
  round <- list()
  
  for (i in seq_along(list)) {
    
    mh <- list[[i]][["MedlineCitation"]][["MeshHeadingList"]]
    
    currentPmid <- list[[i]][['MedlineCitation']][['PMID']][['text']] %>%
      as.numeric()
    
    if (is.null(mh)) {
      
      mhs <- data.frame(
        pmid = currentPmid, muid = NA, mjr = NA, quid = NA, qMjr = NA
      )
      
    } else {
      
      mhs <- data.frame(
        pmid = currentPmid,
        muid = fnM('UI'),
        mjr = fnM('MajorTopicYN'),
        quid = fnQ('UI'),
        qMjr = fnQ('MajorTopicYN'),
        row.names = NULL
      ) %>% separate_rows(all_of(c('quid', 'qMjr')), sep = ";")
      
    }
    
    round[[i]] <- mhs
    
  }
  
  results[[j]] <- round %>% bind_rows
  
  j <- j + 1
  
  if (max == length (pmids)) break
  min <- min + 50
  max <- max + 50
  if (max>length(pmids)) max <- length(pmids)
  
  sleep(0.1)
}

saveRDS(results, file = 'rentrezResults.rds')

resultsDf <- results %>% bind_rows() %>% 
  mutate(
    qMjr = ifelse(is.na(qMjr), 'N', qMjr),
    mjr = ifelse(mjr == 'Y' | qMjr == 'Y', T, F),
    w2 = ifelse(mjr == T, 2, 1),
    w3 = ifelse(mjr == T, 3, 1)
  ) %>% 
  filter(
    !is.na(muid) & #Remove commentaries, corrigendums, etc.
    !(muid %in% c("D005260", "D008297")) #Remove check tags
  ) %>% 
  select(-qMjr)
saveRDS(resultsDf, file = 'mh6PDQ.rds')

pmidsNotRetrieved <- setdiff(pmids, resultsDf$pmid)
save(pmidsNotRetrieved, file = 'pmidsNotRetrieved6.rda')
pmidsNoMH <- resultsDf %>% filter(is.na(muid))
save(pmidsNoMH, file = 'pmidsNoMH6.rda')

mh6 <- resultsDf %>% 
  select(-quid) %>%
  group_by(pmid, muid) %>% 
  mutate( #because we remove P-M duplicates
    mjr = any(mjr),
    w2 = max(w2),
    w3 = max(w3)
  ) %>% 
  distinct(pmid, muid, .keep_all = T)
saveRDS(mh6, file = 'mh6.rds')

nMh6 <- mh6 %>% group_by(pmid) %>% summarise(n = n())
saveRDS(nMh6, file = 'nMh6.rds')

# ***************************************************************** #
#*******************************************************************#
### Information Content                                             #
#*******************************************************************#
# ################################################################# #

mh6 <- readRDS("mh6.rds")
load("tree.rda")

mh6_freq <- mh6 %>% 
  group_by(muid) %>%
  summarize(n = n()) %>%
  right_join(tree, by = 'muid') %>%
  select(muid, n, desc)

desc <- mh6_freq %>%
  filter(!is.na(desc)) %>%
  separate_rows(desc, sep = ';') %>%
  left_join(mh6_freq, by = c('desc' = 'muid')) %>% 
  select(muid, n.y) %>%
  group_by(muid) %>% summarise(ndesc = sum(n.y, na.rm = T))

mh6Ic <- left_join(mh6_freq, desc, by = 'muid') %>%
  select(-desc) %>%
  mutate(
    n     = ifelse(is.na(n), 0, n),
    ndesc = ifelse(is.na(ndesc), 0, ndesc),
    ntot  = n + ndesc,
    ic    = -log(ntot / sum(ntot)) %>% ifelse(is.infinite(.), NA, .)
  )

saveRDS(mh6Ic, file = "mh6Ic.rds")

# ****************************************************************** #
#********************************************************************#
### Distances between terms                                          #
#********************************************************************#
# ****************************************************************** #

# ********************************************* #
#  Step length is 1                             #
# ********************************************* #

load("edgelist.rda")

elClean <- edgelist %>% select(muid, chld)

g1 <- graph_from_data_frame(elClean, directed = F)

dm1 <- distances(g1)

saveRDS(dm1, file = 'dm1.rds', compress = F)

# ********************************************* #
#  Step length is delta-IC                      #
# ********************************************* #

mh6Ic <- readRDS('mh6Ic.rds')

edgelist6T1 <- edgelist %>%
  left_join(mh6Ic,by = c('chld' = 'muid')) %>% #Add IC of children
  select(-ntot, -n, -ndesc) %>% rename(cIc = ic) %>%
  left_join(mh6Ic, by = 'muid') %>% #Add IC of parents
  select(-tn) %>%
  filter(!is.na(cIc)) %>%
  arrange(ic)

#Add IC for categories
cats <- edgelist6T1 %>%
  filter(str_detect(muid, '^[A-Z]$')) %>%
  mutate(
    ndesc = sapply(1:n(), \(x) {
      row <- which(edgelist6T1$muid == .$chld[x]) %>% .[1]
      ntot <- edgelist6T1$ntot[row]
   }),
    ntot = ave(ndesc, muid, FUN = sum),
    ic = -log(ntot/sum(mh6Ic$ntot)),
    n = 0
  ) %>%
  arrange(muid)

edgelist6T2 <- edgelist6T1 %>%
  filter(!str_detect(muid, '^[A-Z]$')) %>% #Remove categories
  rbind(cats) %>%                          #... and then add them with IC
  select(muid, ic, chld, cIc) %>%
  mutate(
    cIc = ifelse( #Add IC to children of ROOT, i.e., categories
      muid == 'ROOT', 
      ic[match(chld, muid)], 
      cIc
    ),
    ic = ifelse(muid == 'ROOT', 0, ic) #Replace ROOT IC with 0
  ) %>%
  arrange(ic)

#Add ROOT, i.e., IC for ROOT and each of 16 children (categories)
tree <- edgelist %>%
  filter(muid == 'ROOT') %>%
  mutate(
    ic = 0,
    cIc = sapply(chld, \(x) {
      row <- which(edgelist6T2$muid == x) %>% .[1]
      edgelist6T2$ic[row]
    })
  ) %>%
  filter(!is.na(cIc)) %>%
  select(muid, ic, chld, cIc)

#Finalise
edgelist6 <- rbind(edgelist6T2, tree) %>% 
  mutate(
    ic = ifelse(is.na(ic), cIc, ic),
    deltaIc = abs(cIc - ic)
  )

save(edgelist6, file = 'edgelist6.rda')

el6Clean <- edgelist6 %>% select(muid, chld, deltaIc)

gIc <- graph_from_data_frame(el6Clean, directed = F)

dmIc <- distances(gIc, weight = E(gIc)$deltaIc) 

saveRDS(dmIc, file = 'dmIc.rds', compress = F)

# ***************************************************************** #
#*******************************************************************#
### Distances between publications                                  #
#*******************************************************************#
# ***************************************************************** #


# ********************************************* #
#  Ready data                                   #
# ********************************************* #

mh <- readRDS('mh6.rds')

#Select relevant PMIDS
rjs <- read.csv('rjs6.txt')
load('rjStats.rda')

pmidsWithMh <- mh %>% pull(pmid) %>% unique()
save(pmidsWithMh, file = 'pmidsWithMh.rda')

topics <- rjStats %>% filter(noRat<0.9) %>% pull(topic)
save(topics, file = 'topics.rda')

pmids <- rjs %>% 
  filter( #only compute for PMIDS:
    pmid %in% unique(mh$pmid) & #... with MH. E.g. 10675423 has RJ but no MH
    topic %in% topics            #... in topics with 10% 'relevant' judgements 
  ) %>% 
  pull(pmid) %>%
  unique() #PMIDS can have RJS for several topics
save(pmids, file = 'pmidsTopics.rda')

pmidList <- list()

for (i in 1:(length(pmids)-1)) {
  sublist <- list(pmids[i], pmids[(i+1):length(pmids)])
  pmidList[[i]] <- sublist
}

save(pmidList, file = 'pmidList.rda')

#Major
mhMjr <- mh %>% 
  filter(pmid %in% pmids) %>% 
  select(-mjr)
saveRDS(mhMjr, file = 'mhMjr.rds')

#Slim
mhSlim <- mh %>% 
  filter(mjr == T & pmid %in% pmids) %>% 
  select(pmid, muid)
saveRDS(mhSlim, file = 'mhSlim.rds')

namesSlim <- unique(mhSlim$muid)

dmSlimIc <- readRDS('dmIc.rds') %>% 
  .[namesSlim, namesSlim]
saveRDS(dmSlimIc, file = 'dmSlimIc.rds')

dmSlim1 <- readRDS('dm1.rds') %>% 
  .[namesSlim, namesSlim]
saveRDS(dmSlim1, file = 'dmSlim1.rds')

#Full
mhFull <- mh %>% 
  select(pmid, muid) %>% 
  filter(pmid %in% pmids)
saveRDS(mhFull, file = 'mhFull.rds')

dmFullIc <- readRDS('dmIc.rds') %>% 
  .[unique(mhFull$muid), unique(mhFull$muid)]
saveRDS(dmFullIc, file = 'dmFullIc.rds')

dmFull1 <- readRDS('dm1.rds') %>%
  .[unique(mhFull$muid), unique(mhFull$muid)]
saveRDS(dmFull1, file = 'dmFull1.rds')

# ********************************************* #
#  Functions                                    #
# ********************************************* #

fnDistS <- function(x) {
  
  focal <- mh %>% filter(pmid == pmidList[[x]][[1]])
  nFocal <- length(unique(focal$muid))
  
  others <- mh %>% filter(pmid %in% pmidList[[x]][[2]])
  
  if (nFocal>1) {
  
    dist <- dm[unique(others$muid), unique(focal$muid)] %>% 
      cbind(., dist = rowMins(.)) %>% 
      data.frame() %>% 
      mutate(muid = rownames(.)) %>% 
      inner_join(others, by = 'muid') %>% 
      group_by(pmid) %>% 
      summarise(
        dnf = sum(dist),
        nNeigh = n(),
        across(.cols = 1:nFocal, .fns = min)
      ) %>% 
      mutate(
        dfn = rowSums(select(., 4:ncol(.))), #Must be in first mutate.
        focal = pmidList[[x]][[1]],
        dist = (dnf + dfn) / (nNeigh + nFocal)
      ) %>% 
      select(pubA = pmid, pubB = focal, dist)
  
  } else if (nFocal == 1) {

    dist <- dm[unique(others$muid), focal$muid] %>% 
      data.frame(dist = .) %>% 
      mutate(muid = rownames(.)) %>% 
      inner_join(others, by = 'muid') %>% 
      group_by(pmid) %>% 
      summarise(
        dnf = sum(dist),
        nNeigh = n(),
        dfn = min(dist)
      ) %>% 
      mutate(
        focal = pmidList[[x]][[1]],
        dist = (dnf + dfn) / (nNeigh + nFocal)
      ) %>% 
      select(pubA = pmid, pubB = focal, dist)
  
  }
  
} 

fnDistW <- function(x) {
  
  focal <- mh %>% filter(pmid == pmidList[[x]][[1]])
  mhF <- unique(focal$muid)
  nFocal <- length(mhF)
  sw2f <- sum(focal$w2)
  sw3f <- sum(focal$w3)
  
  others <- mh %>% filter(pmid %in% pmidList[[x]][[2]])
  
  dist <- dm[unique(others$muid), mhF] %>%
    cbind(., dist = rowMins(.)) %>% 
    data.frame() %>%
    mutate(muid = rownames(.)) %>% 
    inner_join(others, by = 'muid') %>% 
    mutate(distW2 = dist*w2, distW3 = dist*w3) %>% 
    group_by(pmid) %>%
    summarise(
      dnfw1 = sum(dist),
      dnfw2 = sum(distW2),
      dnfw3 = sum(distW3),
      sw2n = sum(w2),
      sw3n = sum(w3),
      nOther = n(),
      across(.cols = 1:nFocal, .fns = min)
    ) 
  
  dist[paste0('TWO_', mhF)] <- lapply(
    mhF,
    \(y) dist[, y] * focal$w2[match(y, focal$muid)]
  )
  
  dist[paste0('THREE_', mhF)] <- lapply(
    mhF,
    \(y) dist[, y] * focal$w3[match(y, focal$muid)]
  )
  
  dist <- dist %>% 
    mutate(
      focal = pmidList[[x]][[1]],
      dfnw1 = rowSums(select(., starts_with('D0'))),
      dfnw2 = rowSums(select(., starts_with('TWO_'))),
      dfnw3 = rowSums(select(., starts_with('THREE_'))),
      distW1 = (dfnw1 + dnfw1) / (nFocal + nOther),
      distW2 = (dfnw2 + dnfw2) / (sw2f + sw2n),
      distW3 = (dfnw3 + dnfw3) / (sw3f + sw3n)
    ) %>%
    select(pubA = pmid, pubB = focal, distW1, distW2, distW3)
  
}


# ********************************************* #
# Delta-IC is step length                       #
# ********************************************* #

rm(list=setdiff(ls(), c("fnDistW", "fnDistS")))

load('pmidList.rda')
load('pmidsTopics.rda')

plan(multisession, workers = 6) 
options(future.globals.maxSize = 800 * 1024^2)

#Slim
dm <- readRDS('dmSlimIc.rds')
mh <- readRDS('mhSlim.rds')

list <- future_lapply(seq_along(pmidList), fnDistS) 

distPubs <- list %>% bind_rows()

saveRDS(distPubs, file = 'distPubs_ic_slim.rds', compress = F)

rm(list, distPubs)

#Full + Major 2 and Major 3
dm <- readRDS('dmFullIc.rds')
mh <- readRDS('mhMjr.rds')

list <- future_lapply(seq_along(pmidList), fnDistW)

distPubs <- list %>% bind_rows()

slim <- readRDS('distPubs_ic_slim.rds')

distPubsIc <- distPubs %>% 
  inner_join(slim, by = c('pubA', 'pubB')) %>% 
  rename(distW1Ic = distW1, distW2Ic = distW2, distW3Ic = distW3,
         distW0Ic = dist)

saveRDS(distPubsIc, file = 'distPubs_ic.rds')

# ********************************************* #
# 1 is step length                              #
# ********************************************* #

rm(list=setdiff(ls(), c("fnDistW", "fnDistS")))

load('pmidList.rda')
load('pmidsTopics.rda')

#Slim
dm <- readRDS('dmSlim1.rds')
mh <- readRDS('mhSlim.rds')

list <- future_lapply(seq_along(pmidList), fnDistS) 

distPubs <- list %>% bind_rows()

saveRDS(distPubs, file = 'distPubs_1_slim.rds', compress = F)

#Full + Major 2 + Major 3
dm <- readRDS('dmFull1.rds')
mh <- readRDS('mhMjr.rds')

list <- future_lapply(seq_along(pmidList), fnDistW)

distPubs <- list %>% bind_rows()

slim <- readRDS('distPubs_1_slim.rds')

distPubs1 <- distPubs %>% 
  inner_join(slim, by = c('pubA', 'pubB')) %>% 
  rename(distW1_1 = distW1, distW2_1 = distW2, distW3_1 = distW3, 
         distW0_1 = dist)

saveRDS(distPubs1, file = 'distPubs_1.rds')

# ********************************************* #
# Co-occurence, Boudreau                        #
# ********************************************* #

mh <- readRDS('mhFull.rds')

#Distance matrix
mat <- sparseMatrix(
  i = match(mh$pmid, unique(mh$pmid)),
  j = match(mh$muid, unique(mh$muid)),
  x = 1,
  dimnames = list(unique(mh$pmid), unique(mh$muid))
) %>% 
  simil(margin = 1, method = 'cosine') %>% 
  as.matrix()

#Distance list
lti <- which(lower.tri(mat), arr.ind = T)

distPubs <- data.frame(
  pubA = as.numeric(rownames(mat)[lti[, 1]]),
  pubB = as.numeric(rownames(mat)[lti[, 2]]),
  distBou = as.numeric(mat[lti])
) %>% 
  mutate( #For joining with similarity distances. Distance is symmetrical
    A = pmax(pubA, pubB),
    B = pmin(pubA, pubB)
  ) %>% 
  select(distBou, pubA = A, pubB = B)

saveRDS(distPubs, file = 'distPubs_boudreau.rds')

# ********************************************* #
# Co-occurence, Ahlgren                         #
# ********************************************* #

load('pmidsTopics.rda')

#Create vector representation
ic <- readRDS("mh6Ic.rds") %>% select(muid, ic)

mh <- readRDS('mh6PDQ.rds') %>% 
  filter(pmid %in% pmids) %>% 
  left_join(ic, by = 'muid') %>% 
  mutate(ic = ifelse(mjr==T, ic*2, ic))

quids <- mh %>% 
  filter(!is.na(quid)) %>% 
  mutate(muid = paste0(muid, ';', quid), ic = 1) %>% 
  select(pmid, muid, ic)

mhQ <- mh %>%
  select(pmid, muid, ic) %>%
  bind_rows(quids)

#Distance matrix
mat <- sparseMatrix(
  i = match(mhQ$pmid, unique(mhQ$pmid)),
  j = match(mhQ$muid, unique(mhQ$muid)),
  x = mhQ$ic,
  dimnames = list(unique(mhQ$pmid), unique(mhQ$muid))
) %>% 
  simil(margin = 1, method = 'cosine') %>% 
  as.matrix()

#Distance list
lti <- which(lower.tri(mat), arr.ind = T)

distPubs <- data.frame(
  pubA = as.numeric(rownames(mat)[lti[, 1]]),
  pubB = as.numeric(rownames(mat)[lti[, 2]]),
  distAhl = as.numeric(mat[lti])
) %>% 
  mutate( #For joining with similarity distances. Distance is symmetrical
    A = pmax(pubA, pubB),
    B = pmin(pubA, pubB)
  ) %>% 
  select(distAhl, pubA = A, pubB = B)

saveRDS(distPubs, file = 'distPubs_ahlgreen.rds')

# ********************************************* #
### Combine                                     #
# ********************************************* #

distPubsAll <- readRDS('distPubs_ic.rds') %>% 
  inner_join(readRDS('distPubs_1.rds'), by = c('pubA', 'pubB')) %>% 
  inner_join(readRDS('distPubs_ahlgreen.rds'), by = c('pubA', 'pubB')) %>% 
  inner_join(readRDS('distPubs_boudreau.rds'), by = c('pubA', 'pubB'))

saveRDS(distPubsAll, file = 'distPubsAll.rds')

load('topics.rda')
distPubs <- readRDS('distPubsAll.rds')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>%
  filter(pmid %in% pmidsWithMh)

fn <- function(x) {
  
  rjsRound <- rjs %>% filter(topic == topics[x])
  
  list <- distPubs %>% 
    inner_join(rjsRound, by = c('pubA' = 'pmid')) %>% 
    inner_join(rjsRound, by = c('pubB' = 'pmid')) %>% 
    mutate(
      rj1 = pmin(rj.x, rj.y),
      rj2 = pmax(rj.x, rj.y)
    ) %>% 
    select(-topic.y, -rj.x, -rj.y, topic = topic.x)
}

list <- lapply(seq_along(topics), fn)

df <- list %>% bind_rows()

saveRDS(df, file = 'bmData.rds')

# ***************************************************************** #
#*******************************************************************#
### Benchmark                                                       #
#*******************************************************************#
# ***************************************************************** #

bmData <- readRDS('bmData.rds')

# ********************************************* #
# Test 1 + Summary stats                        #
# ********************************************* #

df02 <- bmData %>% filter(rj1 == 0 & rj2 == 2)
df22 <- bmData %>% filter(rj1 == 2 & rj2 == 2)

vars <- colnames(bmData) %>% setdiff(c('pubA', 'pubB', 'topic', 'rj1', 'rj2'))

clfD <- future_lapply(
  vars,
  \(x) cliff.delta(df02[[x]], df22[[x]])$estimate
) %>% unlist()

results <- data.frame(
  var = vars,
  mean02 = colMeans(df02[vars]),
  mean22 = colMeans(df22[vars]),
  median02 = apply(df02[vars], 2, median),
  median22 = apply(df22[vars], 2, median),
  sd02 = apply(df02[vars], 2, sd),
  sd22 = apply(df22[vars], 2, sd),
  clfD
) %>% mutate(across(where(is.numeric), ~round(., digits = 3)))

saveRDS(results, file = 'results.RDS')

# ********************************************* #
# Test 2                                        #
# ********************************************* #

load('topics.rda')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>%
  filter(pmid %in% pmidsWithMh)

resultsAll <- list()

for (j in seq_along(topics)) {
  
  print(j)
  
  topicNumber = topics[j]
  
  df <- bmData %>% filter(topic == topicNumber)
  
  pmids <- union(df$pubA, df$pubB) 
  
  rjsTopic <- rjs %>% filter(topic == topicNumber & pmid %in% pmids)
  
  nDecision <- 0
  
  resultsTopic <- data.frame(
    metric = vars,
    tp = 0,
    tn = 0,
    fp = 0,
    fn = 0
  )
  
  for (i in 1:30) {
    
    print(i)
    
    set.seed(i)
    topicPmids <- rjsTopic %>% 
      filter(rj==2) %>% pull(pmid) %>% sample(size = 10)
    
    set.seed(i)
    notTopicPmids <- rjsTopic %>% 
      filter(rj==0) %>% pull(pmid) %>% sample(size = 10)
    
    considerPmids <- rjsTopic %>% 
      filter(!pmid %in% c(topicPmids, notTopicPmids)) %>% pull(pmid)
    
    for (consider in considerPmids) {
      
      resultsTopic$truth <- rjsTopic$rj[rjsTopic$pmid == consider]
      
      if (resultsTopic$truth[1]==1) next
      
      test <- df %>% filter(pubA == consider | pubB == consider)
      
      topic <- test %>% filter(pubA %in% topicPmids | pubB %in% topicPmids)
      
      notTopic <- test %>% filter(pubA %in% notTopicPmids | pubB %in% notTopicPmids)
      
      #If distance (similarity) to pmids in topic is smallest, then pmid under consideration belongs to topic.
      resultsTopic$decision <- sapply(vars, function(var) {
        ifelse(
          str_detect(var, 'W'),
          ifelse(min(notTopic[[var]]) > min(topic[[var]]), 2, 0), 
          ifelse(min(notTopic[[var]]) < min(topic[[var]]), 2, 0)
        ) #if the 'not topics' are further (less similar) than the topics, decision is '2 - relevant'
      })
      
      resultsTopic <- resultsTopic %>% 
        mutate(
          tp = ifelse(truth==2 & decision == 2, tp+1, tp),
          tn = ifelse(truth==0 & decision == 0, tn+1, tn),
          fp = ifelse(truth==0 & decision == 2, fp+1, fp),
          fn = ifelse(truth==2 & decision == 0, fn+1, fn)
        )
      
      nDecision <- nDecision + 1
      
    }
    
  }
  
  resultsTopic <- resultsTopic %>% 
    mutate(
      precision   = round(tp / (tp + fp), 3),
      recall      = round(tp / (tp + fn), 3),
      specificity = round(tn / (tn + fp), 3),
      sensitivity = round(tp / (tp + fn), 3),
      topic       = topicNumber
    ) %>% select(-truth, -decision)
  
  resultsAll[[j]] <- resultsTopic
  
}

information <- resultsAll %>% bind_rows() %>% 
  filter(!str_detect('sim', metric)) %>% 
  mutate(n = tp+tn+fp+fn) %>% 
  group_by(metric) %>% 
  summarize(
    tp = sum(tp),
    tn = sum(tn),
    fp = sum(fp),
    fn = sum(fn)
  ) %>% 
  mutate(
    n = tp+tn+fp+fn,
    precision = round(tp / (tp+fp), 3),
    recall = round(tp / (tp + fn), 3),
    phi = ((tp*tn)-(fp*fn))/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
  )

write.xlsx(information, "informationShort.xlsx")

saveRDS(information, file = 'information.rds')


# ***************************************************************** #
#*******************************************************************#
### Visualise                                                       #
#*******************************************************************#
# ***************************************************************** #

bmData <- readRDS('bmData.rds')
df02 <- bmData %>% filter(rj1 == 0 & rj2 == 2)
df22 <- bmData %>% filter(rj1 == 2 & rj2 == 2)

fnLong <- function(vars) {
  longNr <- df02 %>%
    select(all_of(vars)) %>% 
    pivot_longer(cols = all_of(vars), names_to = "var", values_to = "val") %>% 
    mutate(type = 'nr')
  
  longRr <- df22 %>%
    select(all_of(vars)) %>% 
    pivot_longer(cols = all_of(vars), names_to = "var", values_to = "val") %>% 
    mutate(type = 'rr')
  
  return(rbind(longNr, longRr))
}

vars <- c('distBou', 'distAhl')
long <- fnLong(vars)
lines <- c("solid", "dashed")

pcooc <- ggplot(long, aes(x = val, color = var, linetype = type)) +
  geom_density(size = 0.5) +
  scale_color_manual(values = viridis_pal()(length(vars)), labels = c('Ahlgreen', 'Boudreau')) +
  scale_linetype_manual(values = lines) +
  labs(linetype = "Type", color = "Method", x = "Similarity", y = "Density") +
  coord_cartesian(xlim = c(0, 0.5)) +
  theme_minimal() +
  theme(legend.position = "top")

#∆IC
vars <- c('distW0Ic', 'distW1Ic', 'distW2Ic', 'distW3Ic')
long <- fnLong(vars)

pIc <- ggplot(long, aes(x = val, color = var, linetype = type)) +
  geom_density(size = 0.5) +
  scale_color_manual(values = viridis_pal()(length(vars)), labels = c("Minor 0", "Major 1", "Major 2", "Major 3")) +
  scale_linetype_manual(values = lines) +
  labs(linetype = 'Type', color = "Weight", x = "Distance (∆IC)", y = "Density") +
  coord_cartesian(xlim = c(0, 17.5), ylim = c(0, 0.45)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 17.5)) +
  theme_minimal() +
  theme(legend.position = "bottom")

#1
vars <- c('distW0_1', 'distW1_1', 'distW2_1', 'distW3_1')
long <- fnLong(vars)

p1 <- ggplot(long, aes(x = val, color = var, linetype = type)) +
  geom_density(size = 0.5) +
  scale_color_manual(values = viridis_pal()(length(vars)), guide = F) +
  scale_linetype_manual(values = lines, guide = F) +
  labs(x = "Distance (1)", y = "Density") +
  coord_cartesian(xlim = c(0, 17.5), ylim = c(0, 0.45)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 17.5)) +
  theme_minimal() 

#Combine
pcooc / pIc / p1

ggsave(
  'plots/combined_v2.png',
  width = 7,
  height = 6.5
)
