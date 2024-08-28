
# This script can be used to replicate the results of the paper 'Measuring
# publication relatedness using controlled vocabularies', which will be
# presented at the 28th International Conference on Science, Technology and
# Innovation Indicators, 2024. 

# Access via: https://doi.org/10.48550/arXiv.2408.15004
 
# Author: Emil Dolmer Alnor, Aarhus University, ea@ps.au.dk

# 0. Packages ####
library(httr)
library(jsonlite)
library(rentrez)
library(xml2)
library(XML)
library(igraph)

library(tidyr)
library(stringr)
library(dplyr)
library(data.table)
library(future.apply)
library(Matrix)
library(matrixStats)

library(effsize)
library(proxyC)

library(patchwork)
library(ggplot2)
library(viridis)
library(plotly)

library(openxlsx)
library(readxl)

# ***************************************************************** #
#*******************************************************************#
# 1. MeSH Hiearchy                                               ####
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

fnDesc <- function(x) {
  
  index <- which(str_detect(tree$tn, x))
  
  nodes <- paste(
    tree$muid[index],
    collapse = ";"
  ) %>%
    ifelse(. == "", NA, .)
}

plan(multisession, workers = 6) 
tree$desc <- future_sapply(tree$pattern, fnDesc)

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
## 1.1 Edgelist                              ####
# ********************************************* #

load("tree.rda")

#We start by creating a data.frame with one row for each node representation of the MeSH-terms.
nodes <- tree %>% separate_rows(tn, sep = ';') %>% select(muid, tn)

#Next we modify this data.frame to show the tree number of the parent for each
#node. Apart from highest lvl MeSH-terms (just below the categories) the parent
#is their tree number with the last digits and '.' removed
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

#We now need to add that the parents of the highest level MeSH-terms are the
#categories
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
# 2. TREC Genomics 06'                                           ####
#*******************************************************************#
# ***************************************************************** #

#topics
download.file(
  "https://dmice.ohsu.edu/trec-gen/data/2006/topics/2006topics.xls",
  destfile = "2006topics.xls",
  mode = "wb"
)

temp <- read_excel("2006topics.xls")

fnTop <- function(type, first, last) {
  
  mat <- temp %>%
    slice(first:last) %>%
    mutate(type = type)
  
  colnames(mat) <- c('nid', 'id', 'gene', 'facet2', 'need', 'type')
  
  return(mat)
}

topics6 <- rbind(
  fnTop("disease", 3, 9),
  fnTop("process", 14, 20),
  fnTop("function", 25, 31),
  fnTop("impact", 36, 42)
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
# 3. Publication MeSH data                                       ####
#*******************************************************************#
# ***************************************************************** #

# An API-key for Entrez is required in this section. The API-key is used to
# speed up the download. To replicate the results without registering for an
# API-key, remove the line where 'ak' is defined and remove', api_key = ak' from
# the 'entrez_fetch'.

# ********************************************* #
##  3.1 Functions                            ####
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
## 3.2 Data and loop                         ####
# ********************************************* #

results <- list()
min <- 1
max <- 50
j <- 1
pmids <- read.csv("pmids6_clean.txt") %>% unlist()
ak <- DEFINE THE API_KEY HERE

while (T == T) {
  
  progress(min, 1:length(pmids))
  start_time <- Sys.time()
  
  pmid <- pmids[min:max]
  
  xml <- entrez_fetch(
    db = "pubmed", id = pmid, rettype = "xml", parsed = T, api_key = ak) 
  
  list <- xmlToList(xml) 
  
  # Error handling: Sometimes the entrez_fetch does not work in the first
  # attempt. In the runs of this script, it has always suceed within 3 attempts.
  # If the maximum attempts is reached, manual error inspection is needed.
  
  attempt <- 2
  
  while(length(list)==1) {
    
    Sys.sleep(1)
    
    cat(paste0('Attempt: ', attempt, '\n'))
    
    xml <- entrez_fetch(
      db = "pubmed", id = pmid, rettype = "xml", parsed = T, api_key = ak) 
    
    list <- xmlToList(xml)
    
    attempt <- attempt+1
    
    if (attempt==10) {
      stop('Max attempts reached')
    }
    
  }
  
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
    mjr = ifelse(mjr == 'Y' | qMjr == 'Y', T, F)
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
  group_by(pmid, muid) %>% #because we remove P-M duplicates
  mutate(mjr = any(mjr)) %>% 
  distinct(pmid, muid, .keep_all = T)
saveRDS(mh6, file = 'mh6.rds')

nMh6 <- mh6 %>% group_by(pmid) %>% summarise(n = n())
saveRDS(nMh6, file = 'nMh6.rds')

# ********************************************* #
## 3.3 Prepare data                          ####
# ********************************************* #

#Select relevant PMIDS
mh <- readRDS('mh6.rds')
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

pmidsDf <- pmids %>% data.frame() %>% rename(., pmid = `.`)
save(pmidsDf, file = 'pmidsDf.rda')

# MeSH  

#Major
mhMjr <- mh %>% semi_join(pmidsDf, by = 'pmid') %>% ungroup()

lapply(2:10, function(x) {
  mhMjr[[paste0('w', x)]] <<- ifelse(mhMjr$mjr, mhMjr$mjr*x, 1)
}
)

mhMjr <- mhMjr %>% mutate(across(where(is.numeric), as.integer))

saveRDS(mhMjr, file = 'mhMjr.rds')

namesFull <- unique(mhMjr$muid)

#Slim
mhSlim <- mhMjr %>% 
  filter(mjr == T) %>% 
  select(pmid, muid)
saveRDS(mhSlim, file = 'mhSlim.rds')

namesSlim <- unique(mhSlim$muid)

save(namesSlim, namesFull, file = 'names.rda')

# ********************************************* #
## 3.4 Information content                   ####
# ********************************************* #

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
# 4. Distances between terms                                      ####
#********************************************************************#
# ****************************************************************** #

# ********************************************* #
## 4.1 Step length is 1                      ####
# ********************************************* #

load("edgelist.rda")

elClean <- edgelist %>% select(muid, chld)

g1 <- graph_from_data_frame(elClean, directed = F)

dm1 <- distances(g1)

saveRDS(dm1, file = 'dm1.rds', compress = F)

rm(g1, dm1, elClean)

# ********************************************* #
## 4.2 Step length is delta-IC               ####
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

rm(el6Clean, tree, edgelist6, gIc, dmIc, edgelist6T2, edgelist6T1, cats, mh6Ic)

# ********************************************* #
## 4.3 Zhu et al.                            ####
# ********************************************* #

### 4.3.1 Common data ####### 

load('tree.rda')

freq <- readRDS(file = 'mhMjr.rds') %>% 
  group_by(muid) %>% summarise(n = n()) %>% select(muid)

load('edgelist6.rda')
cats <- edgelist6 %>% 
  filter(nchar(muid) == 1) %>% 
  distinct(muid, ic) %>% 
  mutate(tn = muid)

tns <- tree %>% 
  select(muid, tn) %>%
  separate_rows(tn, sep = ';') %>% 
  inner_join(readRDS('mh6Ic.rds'), by = 'muid') %>% 
  select(muid, tn, ic) %>% 
  bind_rows(cats)

icm <- tns %>% distinct(muid, .keep_all = T) %>% bind_rows(cats)

save(tns, icm, file = 'tns_icm.rda')

df <- tns %>% 
  inner_join(freq, by = 'muid') %>%
  mutate(tn = str_extract_all(tn, "[A-Z]|\\d+")) %>% 
  select(-ic) %>% 
  group_by(muid) %>% 
  mutate(n_tn = n()) %>% 
  ungroup() 

rm(tree, freq, tns, icm, cats, edgelist6)

fnLca <- function(a, b) {
  if (a[1] != b[1]) return(0) #allows early terminate, since T in most cases
  
  len <- min(length(a), length(b))
  
  i <- 1
  
  while (i <= len && a[i] == b[i]) {
    i <- i + 1
  }
  
  return(i-1)
}

### Common dataframe ********

df2 <- bind_rows(
  df %>% inner_join(df, join_by(muid < muid)), 
  df %>% inner_join(df, join_by(muid > muid))
) %>% 
  mutate(dep_cca = mapply(fnLca, tn.x, tn.y))

fnCa <- function(tn, dep) {
  if (is.na(dep)) return(NA)
  if (dep==1) return(tn[1])
  paste0(tn[1], paste(tn[2:dep], collapse = '.'))
}

df_noca <- df2 %>% filter(dep_cca == 0)
df_ca <- df2 %>% filter(dep_cca > 0)
rm(df2, df)

df_ca <- df_ca %>% 
  mutate(
    dep_tnx = vapply(tn.x, length, numeric(1)),
    dep_tny = vapply(tn.y, length, numeric(1)),
    ca = mapply(fnCa, tn.x, dep_cca)
  )

print(1) #Because computations are time-demanding, we split the process up, and
#monitor progress by printing integers

df_ca <- df_ca %>% 
  bind_rows(df_noca) %>% 
  select(-tn.y) %>% 
  mutate(tn.x = vapply(
    tn.x, function(x) paste0(x, collapse = '.'), character(1)) %>%
      factor() %>% as.numeric()
  )

rm(df_noca)
print(2)

df_ca <- df_ca %>% 
  group_by(tn.x, muid.y) %>% 
  filter(dep_cca == max(dep_cca)) %>% 
  filter(row_number() ==1) %>% 
  mutate(
    muid.x2 = pmin(muid.x, muid.y),
    muid.y2 = pmax(muid.x, muid.y),
    n = n_tn.x + n_tn.y,
  ) %>% 
  select(muid.x = muid.x2, tn.x, dep_tnx, muid.y = muid.y2, dep_tny, dep_cca, 
         ca, n)

df3 <- df_ca
rm(df_ca, fnLca, fnCa)

print(3)

saveRDS(df3, file = 'zhu_df3.rds')

### 4.3.2 Path length ########

path <- df3 %>% 
  select(-ca) %>% 
  mutate(
    wp = (2*dep_cca) / (dep_tnx + dep_tny),
    wp_dist = dep_tnx + dep_tny - 2*dep_cca,
    wp1 = exp(-wp_dist/1),
    wp2 = exp(-wp_dist/2),
    wp3 = exp(-wp_dist/3),
    wp4 = exp(-wp_dist/4),
    wp5 = exp(-wp_dist/5)
  ) %>% 
  group_by(muid.x, muid.y) %>% 
  summarise(across(c(wp, wp1, wp2, wp3, wp4, wp5), ~ sum(.x) / first(n)))  %>% 
  mutate(across(c(wp, wp1, wp2, wp3, wp4, wp5), ~ ifelse(is.na(.), 0, .)))

saveRDS(path, file = 'zhu_path.rds')

#Matrix
mhs <- sort(unique(c(path$muid.x, path$muid.y)))

temp_mat <- matrix(
  NA, 
  nrow = length(mhs), 
  ncol = length(mhs), 
  dimnames = list(mhs, mhs)
)

diag(temp_mat) <- 1

load('names.rda')

fnMat <- function(var, dm) {
  print(var)
  
  mat <- temp_mat
  
  mat[as.matrix(rev(dm[, c('muid.x', 'muid.y')]))] <- dm[[var]]
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  
  mat <- mat[namesFull, namesFull]
  saveRDS(mat, file = paste0('dm_f_', var, '.rds'))
  
  mat <- mat[namesSlim, namesSlim]
  saveRDS(mat, file = paste0('dm_s_', var, '.rds'))
}

for (v in c('wp', 'wp1', 'wp2', 'wp3', 'wp4', 'wp5')) fnMat(v, path)

rm(path)

### 4.3.3 IC #####


load('tns_icm.rda')

#Add ic to tns
ic <- df3 %>% 
  select(-dep_tnx, -dep_tny) %>% 
  filter(dep_cca > 0) %>% #only do computations where its needed
  left_join(tns, by = c('ca' = 'tn')) %>% 
  rename(ca_muid = muid, ca_ic = ic) %>% 
  left_join(icm, by = c('muid.x' = 'muid')) %>% rename(ic.x = ic) %>% 
  left_join(icm, by = c('muid.y' = 'muid')) %>% rename(ic.y = ic) %>%
  select(-ca, -ca_muid, -tn.x.x, -tn.y) %>% 
  bind_rows(df3 %>%
              filter(dep_cca == 0) %>%
              select(-dep_tnx, -dep_tny, -ca)
  ) %>% 
  mutate(
    lin = ifelse(dep_cca == 0, 0, (2*ca_ic) / (ic.x + ic.y)),
    djc = ifelse(dep_cca == 0, 9999, ic.x + ic.y - 2*ca_ic),
    jc1 = ifelse(dep_cca == 0, 0, exp(-djc/1)),
    jc2 = ifelse(dep_cca == 0, 0, exp(-djc/2)),
    jc3 = ifelse(dep_cca == 0, 0, exp(-djc/3)),
    jc4 = ifelse(dep_cca == 0, 0, exp(-djc/4)),
    jc5 = ifelse(dep_cca == 0, 0, exp(-djc/5)),
  ) %>%
  group_by(muid.x, muid.y) %>% 
  summarise(across(c(lin, jc1, jc2, jc3, jc4, jc5), ~ sum(.x) / first(n))) 

rm(icm, tns, df3)

saveRDS(ic, file = 'zhu_ic.rds')

#Matrix
for (v in c('lin', 'jc1', 'jc2', 'jc3', 'jc4', 'jc5')) fnMat(v, ic)

rm(temp_mat, namesSlim, namesFull, mhs, fnMat, ic, v)

# ***************************************************************** #
#*******************************************************************#
# 5. Distances between publications                             #####
#*******************************************************************#
# ***************************************************************** #

# ********************************************* #
## 5.1 Distance matrices                    #####
# ********************************************* #

fnDm <- function(dist, names, file) {
  dm <- dist[names, names]
  dm <- -dm # '-' allows to take rowMaxs for dist
  
  saveRDS(dm, file = paste0(file, '.rds'))
  
  dm <- -dm
  
  for (i in 1:5) {
    cat(file, i)
    
    sim <- (exp(-dm / i)) 
    saveRDS(sim, file = paste0(file, '_sim', i, '.rds'))
  }
}

ic <- readRDS('dmIc.rds')
dm1 <- readRDS('dm1.rds') 
load('names.rda')

fnDm(dist = ic, names = namesFull, file = 'dm_f_ic') 
fnDm(dist = ic, names = namesSlim, file  = 'dm_s_ic')
fnDm(dist = dm1, names = namesFull, file = 'dm_f_1')
fnDm(dist = dm1, names = namesSlim, file = 'dm_s_1')

rm(fnDm, ic, dm1, namesFull, namesSlim)

# ********************************************* #
## 5.2 Functions                            #####
# ********************************************* #

fnSimW <- function(x) {
  
  #Values
  focal <- subset(mh, pmid == pmidsDf$pmid[x])
  mhF <- unique(focal$muid)
  nFocal <- length(mhF)
  sumMjr <- sum(focal$mjr)
  for (i in 2:10) assign(paste0('sw', i, 'f'), sum(focal[[paste0('w', i)]]))
  
  others <- mh %>% semi_join(
    pmidsDf[x+1:nrow(pmidsDf), , drop = F],
    by = 'pmid'
  )
  
  #Get closest mesh in focal
  dist <- dm[unique(others$muid), mhF] %>%
    cbind(., dist = rowMaxs(.)) %>% 
    as.data.table(keep.rownames = 'muid') %>% 
    inner_join(others, by = 'muid' )
  
  for (i in 2:10){
    dist[[paste0('dw', i)]] <- dist[['dist']] * dist[[paste0('w', i)]]
  } 
  
  dist <- dist[
    , c(
      .(dnfw1 = sum(dist), 
      dnfw2 = sum(dw2), dnfw3 = sum(dw3), dnfw4 = sum(dw4), 
      dnfw5 = sum(dw5), dnfw6 = sum(dw6), dnfw7 = sum(dw7), 
      dnfw8 = sum(dw8), dnfw9 = sum(dw9), dnfw10 = sum(dw10),
      sw2n = sum(w2), sw3n = sum(w3), sw4n = sum(w4), sw5n = sum(w5), 
      sw6n = sum(w6), sw7n = sum(w7), sw8n = sum(w8), sw9n = sum(w9), 
      sw10n = sum(w10),
      nOther = .N),
      lapply(.SD, max) #Max because dist is negative
    ), 
    .SDcols = mhF, by = pmid
  ]
  
  #Sum closest in neighbor
  for (i in 2:10) {
    dist[, paste0('W', i, mhF) := lapply(
      mhF,
      \(y) .SD[[y]] * focal[[paste0('w', i)]][match(y, focal$muid)]
    ), .SDcols = mhF]
  }
  
  #Non-weighted distance
  dist[, focal := pmidsDf$pmid[x]]
  dist[
    , distW1 := (rowSums(.SD, na.rm = TRUE) + dnfw1) / (nFocal + nOther),
    .SDcols = mhF
  ]
  
  #Weighted distances
  for (i in 2:10) {
    
    dfnw <- rowSums(dist %>% select(starts_with(paste0('W', i))))
    dnfw <- dist[[paste0('dnfw', i)]]
    swf <- get(paste0('sw', i, 'f'))
    swn <- dist[[paste0('sw', i, 'n')]]
    
    dist[[paste0('distW', i)]] <- (dfnw + dnfw) / (swf + swn)
    
  }
  
  #Clean
  dist <- dist %>% select(pubA = pmid, pubB = focal, distW1:distW10)
  
}

fnSimS <- function(x) {
  
  #Values
  focal <- subset(mh, pmid == pmidsDf$pmid[x])
  mhF <- unique(focal$muid)
  nFocal <- length(mhF)
  
  others <- mh %>% semi_join(
    pmidsDf[x+1:nrow(pmidsDf), , drop = F],
    by = 'pmid'
  )
  
  if (nFocal>1) {
    
    dist <- dm[unique(others$muid), mhF] %>%
      cbind(., dist = rowMaxs(.)) %>% 
      as.data.table(keep.rownames = 'muid', key = c('pmid', 'muid')) %>% 
      inner_join(others, by = 'muid')
    
    dist <- dist[
      , c(
        .(dnf = sum(dist),
          nNeigh = .N),
        lapply(.SD, max) #Max because dist is negative
      ), 
      .SDcols = mhF, 
      by = pmid
    ]
    
    dist[, dfn   := rowSums(.SD), .SDcols = mhF]
    dist[, focal := pmidsDf$pmid[x]]
    dist[, dist  := (dnf + dfn) / (nNeigh + nFocal)]
    
    dist <- dist[, .(pubA = pmid, pubB = focal, dist)]
    
  } else if (nFocal == 1) {
    
    dist <- dm[unique(others$muid), mhF] %>% 
      as.data.table(keep.rownames = T, key = c('rn', '.')) 
    names(dist) <- c('muid', 'dist')
    dist <- dist %>% inner_join(others, by = 'muid')
    
    dist <- dist[
      , c(.(
        dnf = sum(dist),
        nNeigh = .N,
        dfn = max(dist) #Max because dist is negative
      )), 
      by = pmid
    ]
    
    dist[, focal := pmidsDf$pmid[x]]
    dist[, dist  := (dnf + dfn) / (nNeigh + nFocal)]
    
    dist <- dist[, .(pubA = pmid, pubB = focal, dist)]
    
  }
  
}

fnRel <- function(mat) {
  print(mat)
  dm <<- readRDS(paste0(mat, '.rds'))
  
  prefix <- str_remove(mat, 'dm_')
  
  if (str_detect(mat, '_s_')) {
    mh <<- readRDS('mhSlim.rds')
    
    list <- future_lapply(1:(nrow(pmidsDf)-1), fnSimS)
    relPubs <- list %>% bind_rows()
    
    names(relPubs)[3] <- prefix
    
  } else if (str_detect(mat, '_f_')) {
    mh <<- readRDS('mhMjr.rds')
    
    list <- future_lapply(1:(nrow(pmidsDf)-1), fnSimW) 
    relPubs <- list %>% bind_rows()
    
    names(relPubs) <- str_replace(
      names(relPubs), 
      '^distW', 
      paste0(prefix, '_w')
    )
  } 
  
  saveRDS(relPubs, file = paste0('rp_', prefix, '.rds'))
}

# ********************************************* #
## 5.3 Graph + Zhu                          #####
# ********************************************* #

load('pmidsDf.rda')

plan(multisession, workers = 22) 
options(future.globals.maxSize = 800 * 1024^2)

matrices <- paste0(
  'dm_', c('f_ic', 's_ic', 'f_1', 's_1'), 
  rep(c('', paste0('_sim', 1:5)), each = 4)
) %>% c(
  paste0(
    c('dm_f_', 'dm_s_'),
    rep(c(paste0('wp', 1:5), paste0('jc', 1:5), 'wp', 'lin'), each = 2)
  )
)

for (matrix in matrices) fnRel(matrix)

rm(dm, fnSimS, fnSimW, fnRel, matrix, pmidsDf)

# ********************************************* #
## 5.4 Co-occurence, Boudreau              ######
# ********************************************* #

mh <- readRDS('mhMjr.rds')

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
## 5.5 Co-occurence, Ahlgren                 ####
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

rm(ic, mh, quids, mhQ, mat, lti, distPubs, pmids)

# ********************************************* #
## 5.6. Combine                             #####
# ********************************************* #

files <- str_replace(matrices, 'dm', 'rp') %>% 
  paste0(., '.rds') %>% c(., 'distPubs_boudreau.rds')

relPubs <- readRDS('distPubs_ahlgreen.rds')

for (file in files) {
  print(file)
  relPubs <- relPubs %>% inner_join(readRDS(file), by = c('pubA', 'pubB'))
}

saveRDS(relPubs, file = 'relPubsAll.rds')

load('topics.rda')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>% filter(pmid %in% pmidsWithMh)

fnRjs <- function(x) {
  
  rjsRound <- rjs %>% filter(topic == topics[x])
  
  list <- relPubs %>% 
    inner_join(rjsRound, by = c('pubA' = 'pmid')) %>% 
    inner_join(rjsRound, by = c('pubB' = 'pmid')) %>% 
    mutate(
      rj1 = pmin(rj.x, rj.y),
      rj2 = pmax(rj.x, rj.y)
    ) %>% 
    select(-topic.y, -rj.x, -rj.y, topic = topic.x)
}

bmData <- lapply(seq_along(topics), fnRjs) %>% bind_rows()

saveRDS(bmData, file = 'bmData.rds')

rm(relPubs, fnRjs, files, rjs, pmidsWithMh, topics)

# ***************************************************************** #
#*******************************************************************#
# 6. Benchmark                                                   ####
#*******************************************************************#
# ***************************************************************** #

# ********************************************* #
## 6.1 Test 1 + Summary stats                ####
# ********************************************* #

df02 <- bmData %>% filter(rj1 == 0 & rj2 == 2)
df22 <- bmData %>% filter(rj1 == 2 & rj2 == 2)

vars <- colnames(bmData) %>% setdiff(c('pubA', 'pubB', 'topic', 'rj1', 'rj2'))

clfD <- future_lapply(
  vars,
  \(x) cliff.delta(df22[[x]], df02[[x]])$estimate
) %>% unlist()

test1 <- data.frame(
  var = vars,
  mean02 = colMeans(df02[vars]),
  mean22 = colMeans(df22[vars]),
  median02 = apply(df02[vars], 2, median),
  median22 = apply(df22[vars], 2, median),
  sd02 = apply(df02[vars], 2, sd),
  sd22 = apply(df22[vars], 2, sd),
  clfD
) %>% mutate(across(where(is.numeric), ~round(., digits = 3)))

save(test1, file = 'test1.rda')

rm(df02, df22, clfD, test1)

# ********************************************* #
## 6.3 Test 2                               #####
# ********************************************* #

load('topics.rda')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>%
  filter(pmid %in% pmidsWithMh & rj %in% c(0, 2))
vars <- colnames(bmData) %>% setdiff(c('pubA', 'pubB', 'topic', 'rj1', 'rj2'))

resultsAll <- list()

for (j in seq_along(topics)) {
  
  cat('topic', j, '\n')
  
  topicNumber <- topics[j]
  
  df <- bmData %>% filter(topic == topicNumber)
  
  pmids <- union(df$pubA, df$pubB) 
  
  rjsTopic <- rjs %>% filter(topic == topicNumber & pmid %in% pmids)
  topicPmids <- rjsTopic %>% filter(rj==2) %>% pull(pmid)
  notTopicPmids <- rjsTopic %>% filter(rj==0) %>% pull(pmid)
  
  resultsRounds <- list()
  
  for (i in 1:30) {
    
    print(i)
    
    set.seed(i)
    topicSample <- topicPmids %>% sample(size = 10)
    
    set.seed(i)
    notTopicSample <- notTopicPmids %>% sample(size = 10)
    
    considerPmids <- rjsTopic$pmid[
      !rjsTopic$pmid %in% c(topicSample, notTopicSample)
    ] 
    
    resultsRounds[[i]] <- lapply(considerPmids, function(consider) {
      
      test <- subset(df, pubA == consider | pubB == consider)
      topic <- subset(test, pubA %in% topicSample | pubB %in% topicSample)
      notTopic <- subset(test, pubA %in% notTopicSample | pubB %in% notTopicSample
      ) 
      
      result <- cbind(
        var = vars,
        truth = rjsTopic$rj[rjsTopic$pmid == consider],
        decision = sapply(vars, function(var) {
          ifelse(max(notTopic[[var]]) < max(topic[[var]]), 2, 0)
        })
      )
      
      return(result)
      
    })}
  
  resultsAll[[j]] <- lapply(
    resultsRounds, 
    function(il) do.call(rbind, il)
  ) %>%
    do.call(rbind, .) %>% 
    data.frame() %>% 
    mutate(
      tp = (truth == 2 & decision == 2),
      tn = (truth == 0 & decision == 0),
      fp = (truth == 0 & decision == 2),
      fn = (truth == 2 & decision == 0)
    ) %>% 
    group_by(var) %>% 
    summarise(across(c('tp', 'tn', 'fp', 'fn'), sum))
}

test2 <- resultsAll %>% bind_rows() %>% 
  group_by(var) %>% 
  summarise(across(c('tp', 'tn', 'fp', 'fn'), sum)) %>% 
  mutate(
    n = tp+tn+fp+fn,
    precision = round(tp / (tp+fp), 3),
    recall = round(tp / (tp + fn), 3),
    phi = ((tp*tn)-(fp*fn))/(sqrt(
      (as.numeric(tp)+as.numeric(fp))* #Avoids integer overflow
        (as.numeric(tp)+as.numeric(fn))*
        (as.numeric(tn)+as.numeric(fp))*
        (as.numeric(tn)+as.numeric(fn))
      ))
  )

save(test2, file = 'test2.rda')

rm(topicNumber, df, pmids, rjsTopic, topicPmids, notTopicPmids, topicSample,
   notTopicSample, considerPmids, resultsRounds, resultsTopicList, 
   resultsTopic, resultsAll, test2)


# ***************************************************************** #
#*******************************************************************#
# 7. Visualise                                                   ####
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

#IC
vars <- c('distW0Ic', 'distW1Ic', 'distW2Ic', 'distW3Ic')
long <- fnLong(vars)

pIc <- ggplot(long, aes(x = val, color = var, linetype = type)) +
  geom_density(size = 0.5) +
  scale_color_manual(values = viridis_pal()(length(vars)), labels = c("Minor 0", "Major 1", "Major 2", "Major 3")) +
  scale_linetype_manual(values = lines) +
  labs(linetype = 'Type', color = "Weight", x = "Distance (IC)", y = "Density") +
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

rm(bmData, df02, df22, long, p1, pIc, lines, vars, fnLong)
