
# 0. Packages ####
library(data.table)
library(dplyr)
library(e1071)
library(effsize)
library(epiR)
library(future.apply)
library(ggplot2)
library(htmlwidgets)
library(httr)
library(igraph)
library(jsonlite)
library(matrixStats)
library(Matrix)
library(openxlsx)
library(patchwork)
library(plotly)
library(proxyC)
library(readxl)
library(rentrez)
library(stringr)
library(tidyr)
library(viridis)
library(XML)
library(xml2)

#Function used to monitor progress in heavy computations
fnLoop <- function(iterator, pref = NULL, suf = NULL) {
  if (is.null(pref) & is.null(suf)) {
    cat(iterator, format(Sys.time(), "%H:%M:%S"), '\n')
  } else if (is.null(pref) & !is.null(suf)) {
    cat(iterator, suf, format(Sys.time(), "%H:%M:%S"), '\n')
  } else if (is.null(suf)) {
    cat(pref, iterator, format(Sys.time(), "%H:%M:%S"), '\n')
  } else {
    cat(pref, iterator, suf, format(Sys.time(), "%H:%M:%S"), '\n')
  }
}

# **************************************************************************** #
#******************************************************************************#
# 1. MeSH Hiearchy                                                          ####
#******************************************************************************#
# **************************************************************************** #

# ********************************************* #
## 1.1 Download and parse                    ####
# ********************************************* #

# Make a df with 1 row pr. MeSH and cols showing tree-number(s) and direct and
# indirect descendants

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

#Extract tree numbers
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

#Descendants
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

#Children (direct descendants)
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

rm(tnl, tn, tree, i, j, s, t, lines, mh, muid, nrl, fnDesc)

# ********************************************* #
## 1.2 Edgelist                              ####
# ********************************************* #

load("tree.rda")

# Start by creating a data.frame with one row for each tree-number of each
# MeSH-terms.
nodes <- tree %>% separate_rows(tn, sep = ';') %>% select(muid, tn)

# Next modify this data.frame to show the tree number of the parent for each
# node. Apart from highest lvl MeSH-terms (just below the categories) the parent
# is their tree number with the last digits and '.' removed
children <- nodes %>% 
  mutate(
    tnp = ifelse(
      str_detect(tn, '\\.'), 
      str_remove(tn, '\\.\\d+$'), 
      str_extract(tn, '^.')
    )
  ) %>% select(muid, tnp) %>% rename(chld = muid)

# Now add the children to each of the MeSH-terms.
edgelist <- nodes %>% 
  inner_join(children, by = c('tn' = 'tnp')) %>% 
  distinct(chld, muid, .keep_all = T) #Many parent-child-pairs are parent-child-pairs with more than 1 node.

# We now add that the parents of the highest level MeSH-terms are the
# categories
children <- edgelist %>% 
  filter(!str_detect(tn, '\\.')) %>% #All nodes with '.' are not highest lvl
  distinct(muid, tn) %>%  
  mutate(tn = str_extract(tn, '^.')) %>%  #We take the category letter
  rename(chld = muid) %>% 
  mutate(muid = tn)

edgelist <- bind_rows(edgelist, children)
saveRDS(edgelist, 'edgelist.rds')

rm(children, edgelist, nodes, tree)

# **************************************************************************** #
#******************************************************************************#
# 2. TREC Genomics 06'                                                      ####
#******************************************************************************#
# **************************************************************************** #

#Download relevance judgments and PMIDS of articles.

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
  summarise(rj = max(rj), .groups = 'drop')

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

rm(rjs6, rjStats, t6, temp, topics6, error_pmids, pmids6, pmids6_clean, fnTop)

# **************************************************************************** #
#******************************************************************************#
# 3. Publication MeSH data                                                  ####
#******************************************************************************#
# **************************************************************************** #

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

saveRDS(results, 'rentrezResults.rds')

#Dataframe with one row for each PMID-MeSH-Qualifier
results <- readRDS('rentrezResults.rds')
pmids <- read.csv("pmids6_clean.txt") %>% unlist()

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
saveRDS(resultsDf, 'mh6PDQ.rds')

pmidsNotRetrieved <- setdiff(pmids, resultsDf$pmid)
save(pmidsNotRetrieved, file = 'pmidsNotRetrieved6.rda')

#Dataframe with one row for each PMID-MeSH
mh6 <- resultsDf %>% 
  select(-quid) %>%
  group_by(pmid, muid) %>% #because we remove P-M duplicates
  mutate(mjr = any(mjr)) %>% 
  distinct(pmid, muid, .keep_all = T)
saveRDS(mh6, file = 'mh6.rds')

#Shows number of MeSH pr PMID
nMh6 <- mh6 %>% group_by(pmid) %>% summarise(n = n())
saveRDS(nMh6, file = 'nMh6.rds')

rm(mh6, nMh6, pmids, pmidsNoMH, xml, list, attempt, round, NoMH, results, resultsDf, j, max, min, pmids, pmidsNotRetrieved, fnM, fnQ, progress, sleep)

# ********************************************* #
## 3.3 Prepare data                          ####
# ********************************************* #

#Select relevant PMIDS
mh <- readRDS('mh6.rds')
rjs <- read.csv('rjs6.txt')
load('rjStats.rda')

pmidsWithMh <- mh %>% pull(pmid) %>% unique()
save(pmidsWithMh, file = 'pmidsWithMh.rda')

#Select topics with at least 90% relevant or possibly relevant
topics <- rjStats %>% filter(noRat<0.9) %>% pull(topic)
save(topics, file = 'topics.rda')

#Now select relevant PMIDS
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

# Process the MeSH data  

#Major
mhMjr <- mh %>% semi_join(pmidsDf, by = 'pmid') %>% ungroup()

#Weights to use in Maximum Term Similarities
lapply(2:20, function(x) {
  mhMjr[[paste0('w', x)]] <<- ifelse(mhMjr$mjr, mhMjr$mjr*x, 1)
})

mhMjr <- mhMjr %>% mutate(across(where(is.numeric), as.integer))

saveRDS(mhMjr, file = 'mhMjr.rds')

#Full and slim MeSH-terms, used to subset MeSH-MeSH-distance/similarity matrices
namesFull <- unique(mhMjr$muid)

#Slim
mhSlim <- mhMjr %>% 
  filter(mjr == T) %>% 
  select(pmid, muid)
saveRDS(mhSlim, file = 'mhSlim.rds')

namesSlim <- unique(mhSlim$muid)

save(namesSlim, namesFull, file = 'names.rda')

rm(mh, mhMjr, mhSlim, pmidsDf, rjs, rjStats, namesFull, namesSlim, pmids, pmidsWithMh, topics)

# ********************************************* #
## 3.4 Information content                   ####
# ********************************************* #

mh6 <- readRDS("mh6.rds")
load("tree.rda")

#Frequency of each MeSH term
mh6_freq <- mh6 %>% 
  group_by(muid) %>%
  summarize(n = n()) %>%
  right_join(tree, by = 'muid') %>%
  select(muid, n, desc)

#Frequency of descendants of each MeSH term
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
    ntot  = n + ndesc, #"Total" frequency is frequency + frequency of descendants
    ic    = -log(ntot / sum(ntot)) %>% ifelse(is.infinite(.), NA, .)
  )

saveRDS(mh6Ic, file = "mh6Ic.rds")

rm(desc, mh6, mh6_freq, mh6Ic, tree)

# **************************************************************************** #
#******************************************************************************#
# 4. Term similarity                                                        ####
#******************************************************************************#
# **************************************************************************** #

# ********************************************* #
## 4.1 Unweighted graph                      ####
# ********************************************* #

edgelist <- readRDS('edgelist.rds')

elClean <- edgelist %>% select(muid, chld)

g1 <- graph_from_data_frame(elClean, directed = F)

dm1 <- distances(g1)

saveRDS(dm1, file = 'dm1.rds', compress = F)

rm(g1, dm1, elClean, edgelist)

# ********************************************* #
## 4.2 IC-weighted graph                     ####
# ********************************************* #

el <- readRDS('edgelist.rds') %>% select(muid, chld)
mh_ic <- readRDS('mh6Ic.rds')
sum_ntot <- sum(mh_ic$ntot)

#Add IC for children and parents
muids <- el %>%
  left_join(mh_ic, by = c('chld' = 'muid')) %>% #Children
  select(-n, -ndesc) %>% rename(c_ic = ic, c_ntot = ntot) %>%
  left_join(mh_ic, by = 'muid') #Parent

#Calculate IC for categories
cats <- muids %>%
  filter(str_detect(muid, '^[A-Z]$')) %>% 
  group_by(muid) %>% mutate(ntot = sum(c_ntot)) %>% 
  mutate(
    ic = -log(ntot / sum_ntot),
    delta_ic =  c_ic - ic
  ) %>% 
  filter(ntot != 0) #Remove 'Publication format'

#Remove categories from 'muids' and calculate change in IC. Then combine the 2.
muids_clean <- muids %>% 
  filter(!str_detect(muid, '^[A-Z]$')) %>% 
  mutate(delta_ic = abs(c_ic - ic))

el_icw <- bind_rows(muids_clean, cats) %>% 
  select(muid, chld, delta_ic) %>% 
  filter(!is.na(delta_ic))

saveRDS(el_icw, 'edgelist_icw.rds')

#Now calculate distances
g_ic <- graph_from_data_frame(el_icw, directed = F)

dm_ic <- distances(g_ic, weight = E(g_ic)$delta_ic) 

saveRDS(dm_ic, file = 'dm_ic.rds', compress = F)

rm(el, mh_ic, sum_ntot, muids, cats, muids_clean, el_icw, g_ic, dm_ic)

# ********************************************* #
## 4.3 Distance matrices                    #####
# ********************************************* #

fnDm <- function(dist, names, file) {
  dm <- dist[names, names] #We only need to calculate distances between MeSH, which are in the publications. The number of MeSH depends on whether minor terms are dropped.
  dm <- -dm # '-' allows to take rowMaxs for dist in fnSimW and fnSimS
  
  saveRDS(dm, file = paste0(file, '.rds'))
  
  dm <- -dm
  
  for (i in 1:5) { #Similarity matrices depend on value of λ
    cat(file, i)
    
    sim <- (exp(-dm / i)) 
    saveRDS(sim, file = paste0(file, '_sim', i, '.rds'))
  }
}

ic <- readRDS('dm_ic.rds')
dm1 <- readRDS('dm1.rds') 
load('names.rda')

fnDm(dist = ic,  names = namesFull, file = 'dm_f_ic') 
fnDm(dist = ic,  names = namesSlim, file = 'dm_s_ic')
fnDm(dist = dm1, names = namesFull, file = 'dm_f_1')
fnDm(dist = dm1, names = namesSlim, file = 'dm_s_1')

rm(fnDm, ic, dm1, namesFull, namesSlim)

# **************************************************************************** #
#******************************************************************************#
# 5. Publication relatedness                                               #####
#******************************************************************************#
# **************************************************************************** #

# Function to combine relatedness scores with relevance judgements of 
# publications in a publication pair. Used in sections '5.x.1 Combine'.
fnRjs <- function(x) {
  
  rjsRound <- rjs %>% filter(topic == topics[x])
  
  list <- df %>% 
    inner_join(rjsRound, by = c('pmid_A' = 'pmid')) %>% 
    inner_join(rjsRound, by = c('pmid_B' = 'pmid')) %>% 
    mutate(
      rj1 = pmin(rj.x, rj.y), # We don't use which pmid got which rj, just 
      rj2 = pmax(rj.x, rj.y)  # whether the pair is 2-2 or 0-2
    ) %>% 
    select(-topic.y, -rj.x, -rj.y, topic = topic.x)
}

# ********************************************* #
## 5.1 Term similarity                      #####
# ********************************************* #

### 5.1.1 Functions #############################

fnSimW <- function(x) {
  
  #Values
  focal <- subset(mh, pmid == pmidsDf$pmid[x])
  mhF <- unique(focal$muid)
  nFocal <- length(mhF)
  sumMjr <- sum(focal$mjr)
  for (i in 2:20) assign(paste0('sw', i, 'f'), sum(focal[[paste0('w', i)]]))
  
  others <- mh %>% semi_join(
    pmidsDf[x+1:nrow(pmidsDf), , drop = F],
    by = 'pmid'
  )
  
  # 'dm' shows the distance (similarity) between all MeSH-terms. We start out by
  # subsetting 'dm', so the MeSH-terms of the focal pmid are the columns, and
  # the MeSH terms of the other pmids are in the rows. In 'cbind' we find the
  # nearest MeSH in the focal article for each of the other MeSH-terms, and
  # compute that distance. The next step is to merge the pmids of the other
  # articles. This is a 1:m merge, creating a df in long-format. To be able to
  # use inner_join, we first convert to data.frame, and get the key (muid) to
  # join by.
  dist <- dm[unique(others$muid), mhF] %>%
    cbind(., dist = rowMaxs(.)) %>% 
    as.data.table(keep.rownames = 'muid') %>% 
    inner_join(others, by = 'muid' )
  
  # Apply weights
  for (i in 2:20){
    dist[[paste0('dw', i)]] <- dist[['dist']] * dist[[paste0('w', i)]]
  } 
  
  # Now we have many MeSH-duplicates but unique PMID-MeSH in the rows. In w2-w20 
  # we have the weights for that MeSH-term in that PMID. 
  # To get the distance from the neighbor articles to the focal article we need
  # to get 1) the sum of the weighted distances between the neighbors and the
  # focal article, 2) the sum of the weights (we need them to compute the
  # distance in the end). For the unweighted version, this is simply the sum of
  # distances and the number of mesh terms in the neighbor articles.
  # To get the distance from the focal article to the neighbor articles, we
  # first need to find out what mesh-term in each of the other articles is
  # closest the each of the focal articles MeSH-terms, and compute this
  # distance. 
  dist <- dist[
    , c(
      .(dnfw1 = sum(dist), 
        dnfw2 = sum(dw2), dnfw3 = sum(dw3), dnfw4 = sum(dw4), dnfw5 = sum(dw5), 
        dnfw6 = sum(dw6), dnfw7 = sum(dw7), dnfw8 = sum(dw8), dnfw9 = sum(dw9),
        dnfw10 = sum(dw10), dnfw11 = sum(dw11), dnfw12 = sum(dw12),
        dnfw13 = sum(dw13), dnfw14 = sum(dw14), dnfw15 = sum(dw15), 
        dnfw16 = sum(dw16), dnfw17 = sum(dw17), dnfw18 = sum(dw18), 
        dnfw19 = sum(dw19), dnfw20 = sum(dw20),
        sw2n = sum(w2), sw3n = sum(w3), sw4n = sum(w4), sw5n = sum(w5), 
        sw6n = sum(w6), sw7n = sum(w7), sw8n = sum(w8), sw9n = sum(w9), 
        sw10n = sum(w10), sw11n = sum(w11), sw12n = sum(w12), sw13n = sum(w13),
        sw14n = sum(w14), sw15n = sum(w15), sw16n = sum(w16), sw17n = sum(w17),
        sw18n = sum(w18), sw19n = sum(w19), sw20n = sum(w20),
        nOther = .N),
      lapply(.SD, max) 
    ), 
    .SDcols = mhF, by = pmid
  ]
  
  # Apply weights
  for (i in 2:20) {
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
  
  # In the final step we first compute the sum of the weighted distances from
  # the focal article to the neighbor articles. Then we compute the distance as:
  # distance from neighbor to focal + distance from focal to neighbor divided by
  # the sum of their weights. For the unweighted version, this is simply the sum
  # of distances divided the number of MeSH-terms
  for (i in 2:20) {
    
    dfnw <- rowSums(dist %>% select(starts_with(paste0('W', i, 'D'))))
    dnfw <- dist[[paste0('dnfw', i)]]
    swf <- get(paste0('sw', i, 'f'))
    swn <- dist[[paste0('sw', i, 'n')]]
    
    dist[[paste0('distW', i)]] <- (dfnw + dnfw) / (swf + swn)
    
  }
  
  #Clean up
  dist <- dist %>% select(pubA = pmid, pubB = focal, distW1:distW20)
  
}

#See fnSimW for documentation.
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
        lapply(.SD, max) 
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
        dfn = max(dist) 
      )), 
      by = pmid
    ]
    
    dist[, focal := pmidsDf$pmid[x]]
    dist[, dist  := (dnf + dfn) / (nNeigh + nFocal)]
    
    dist <- dist[, .(pubA = pmid, pubB = focal, dist)]
    
  }
  
}

# Function to read the distance matrix, choose correct MeSH-file, calculate
# reledness, rename relatedness columns, and save file.
fnRel <- function(mat) {

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
  rm(list, relPubs, prefix, dm)
}


### 5.1.2 Calculate #############################

load('pmidsDf.rda')

plan(multisession, workers = 6) 
options(future.globals.maxSize = 800 * 1024^2)

matrices <- paste0(
  'dm_', c('f_ic', 's_ic', 'f_1', 's_1'), 
  rep(c('', paste0('_sim', 1:5)), each = 4)
)

for (matrix in matrices) { fnLoop(matrix)
  fnRel(matrix)
}

rm(dm, fnSimS, fnSimW, fnRel, matrix, pmidsDf)
plan(sequential)

### 5.1.3 Combine ###############################

# Combine into 1 df. This is done slightly different than Salton's and Soft 
# Cosine because I was forced to come up with a more efficient method, because
# I only had my laptop (with low RAM) at time of computation.

#Prepare for benchmark
load('topics.rda')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>% filter(pmid %in% pmidsWithMh)

files <- matrices %>% str_replace('dm', 'rp') %>% paste0(., '.rds')

for (file in files) { fnLoop(file)
  df <- readRDS(file) %>% rename(pmid_A = pubA, pmid_B = pubB)
  bm <- lapply(seq_along(topics), fnRjs) %>% bind_rows() %>% setDT()
  rm(df)
  setkey(bm, pmid_A, pmid_B)
  bmfile <- file %>% str_replace('rp', 'bm')
  saveRDS(bm, bmfile)
}
rm(bm)

bmfiles <- str_replace(files, 'rp', 'bm')
ts_bm <- readRDS(bmfiles[1])
for (file in bmfiles[-1]) { fnLoop(file) 
  ts_bm <- merge(ts_bm, readRDS(file), by = c('pmid_A', 'pmid_B', 'rj1', 'rj2', 'topic'))
}

saveRDS(ts_bm, 'ts_bm.rds')

df02 <- ts_bm %>% filter(rj1 == 0 & rj2 == 2)
saveRDS(df02, file = 'ts_02.rds')
df22 <- ts_bm %>% filter(rj1 == 2 & rj2 == 2)
saveRDS(df22, file = 'ts_22.rds')

rm(df, file, files, matrices, rjs, pmidsWithMh, topics, df02, df22, ts_bm)

# ********************************************* #
## 5.2 Salton's cosine                     ######
# ********************************************* #

load('pmidsTopics.rda')

### 5.2.1 Calculate #############################

# Calculate cosine similarity using sparse matrices of pmid-muid. Then convert
# to data.tables for fast merging.

fnMat2dt <- function(matrix, colname) {
  
  df <- data.table(
    tA = as.integer(rownames(matrix)[lti[, 1]]),
    tB = as.integer(rownames(matrix)[lti[, 2]]),
    sim = as.numeric(matrix[lti])
  )
  df[, 
     `:=`(pmid_A = pmin(tA, tB), 
          pmid_B = pmax(tA, tB))
  ][, c("tA", "tB") := NULL]
  
  setnames(df, 'sim', name) 
  setkey(df, 'pmid_A', 'pmid_B')
  
  return(df)
}

#Without qualifiers
mh <- inner_join(
  readRDS('mhMjr.rds') %>% select(pmid, muid, mjr),
  readRDS('mh6Ic.rds') %>% select(muid, ic),
  by = 'muid'
) %>% filter(pmid %in% pmids)

is <- match(mh$pmid, unique(mh$pmid))
js <- match(mh$muid, unique(mh$muid))
dns <- list(unique(mh$pmid), unique(mh$muid))

for (major_weight in 1:10) { print(major_weight)
  
  name <- paste0('cc_1_w', major_weight)
  
  # 1
  mat <- sparseMatrix(
    i = is,
    j = js,
    x = ifelse(mh$mjr, major_weight, 1),
    dimnames = dns,
  ) %>% 
    simil(margin = 1, method = 'cosine') %>% 
    as.matrix()
  
  if (major_weight == 1) lti <- which(lower.tri(mat), arr.ind = T)
  
  df <- fnMat2dt(matrix = mat, colname = name)
  saveRDS(df, paste0(name, '.rds'))
  
  # IC
  name <- str_replace(name, '_1_', '_ic_')
  
  mat <- sparseMatrix(
    i = is,
    j = js,
    x = ifelse(mh$mjr, major_weight * mh$ic, mh$ic),
    dimnames = dns,
  ) %>% 
    simil(margin = 1, method = 'cosine') %>% 
    as.matrix()
  
  df <- fnMat2dt(matrix = mat, colname = name)
  saveRDS(df, paste0(name, '.rds'))
}

#With qualifiers
quids <- readRDS('mh6PDQ.rds') %>% 
  filter(!is.na(quid) & pmid %in% pmids) %>% 
  mutate(muid = paste0(muid, ';', quid), ic = 1, mjr = F) %>% 
  select(pmid, muid, mjr, ic)

mhq <- bind_rows(quids, mh)

is <- match(mhq$pmid, unique(mhq$pmid))
js <- match(mhq$muid, unique(mhq$muid))
dns <- list(unique(mhq$pmid), unique(mhq$muid))

for (major_weight in 1:10) { cat(major_weight, ' ', format(Sys.time(), "%H:%M:%S"), '\n')
  
  name <- paste0('cc_ic_q_w', major_weight)
  
  mat <- sparseMatrix(
    i = is,
    j = js,
    x = ifelse(mhq$mjr, major_weight * mhq$ic, mhq$ic),
    dimnames = dns
  ) %>% 
    simil(margin = 1, method = 'cosine') %>% 
    as.matrix()
  
  if (major_weight == 1) lti <- which(lower.tri(mat), arr.ind = T)
  
  df <- fnMat2dt(matrix = mat, colname = name)
  saveRDS(df, paste0(name, '.rds'))
}

rm(df, dns, lti, mat, mh, mhq, quids, is, js, major_weight, name, pmids, fnMat2dt)

### 5.2.2 Combine ###############################

#Combine
cc_files <- c(
  paste0('cc_ic_w', 1:10, '.rds'),
  paste0('cc_1_w', 1:10, '.rds'),
  paste0('cc_ic_q_w', 1:10, '.rds')
)

df <- readRDS(cc_files[1])

for (cc_file in cc_files[-1]) { print(cc_file)
  df <- merge(df, readRDS(cc_file), by = c('pmid_A', 'pmid_B'))
}

saveRDS(df, 'cc.rds')

#Add topics and relevance judgement
load('topics.rda')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>% filter(pmid %in% pmidsWithMh)

bm_data <- lapply(seq_along(topics), fnRjs) %>% bind_rows()
saveRDS(bm_data, file = 'cc_bm.rds')

cc_02 <- bm_data %>% filter(rj1 == 0 & rj2 == 2)
saveRDS(cc_02, 'cc_02.rds')

cc_22 <- bm_data %>% filter(rj1 == 2 & rj2 == 2)
saveRDS(cc_22, 'cc_22.rds')

# Compute concordance correlation coefficient to measure agreement between 
# relatedness measures

ccc <- data.frame(icVicq = rep(NA, 10), icV1 = rep(NA, 10), w = 1:10)

# Compare vector with and without qualifiers
ccc$icVicq <- sapply(1:10, function(i) {
  epi.ccc(
    bm_data[[paste0('cc_ic_w', i)]], 
    bm_data[[paste0('cc_ic_q_w', i)]]
  )$rho.c$est
})

# Compare to weighted and non-weighted vectors, which should not agree.
ccc$icV1 <- sapply(1:10, function(i) {
  epi.ccc(
    bm_data[[paste0('cc_ic_w', i)]], 
    bm_data[[paste0('cc_1_w', i)]]
  )$rho.c$est
})

print(print(ccc, digits = 8))

saveRDS(ccc, 'ccc.rds')

rm(bm_data, cc_02, cc_22, df, rjs, cc_file, cc_files, pmidsWithMh, topics, ccc)


# ********************************************* #
## 5.3 Soft cosine                           ####
# ********************************************* #

### 5.3.1 Calculate #############################

# Create vectors representation of pmids depending on 1) weight given to major 
# terms and 2) whether we weight by information content or not. Store in lists
# for fast acess.
mh <- inner_join(
  readRDS('mhMjr.rds') %>% select(pmid, muid, mjr),
  readRDS("mh6Ic.rds") %>% select(muid, ic),
  by = 'muid'
)

for (major_weight in 1:10) {  print(major_weight)
  
  #1
  pm_list <- sparseMatrix(
    i = match(mh$pmid, unique(mh$pmid)),
    j = match(mh$muid, unique(mh$muid)),
    x = ifelse(mh$mjr, major_weight, 1),
    dimnames = list(unique(mh$pmid), unique(mh$muid))
  ) %>% as.matrix() %>% split(., seq_len(nrow(.)))
  
  saveRDS(pm_list, paste0('pm_list_1_', major_weight, '.rds'))
  
  #IC
  pm_list <- sparseMatrix(
    i = match(mh$pmid, unique(mh$pmid)),
    j = match(mh$muid, unique(mh$muid)),
    x = ifelse(mh$mjr, mh$ic*major_weight, mh$ic),
    dimnames = list(unique(mh$pmid), unique(mh$muid))
  ) %>% as.matrix() %>% split(., seq_len(nrow(.)))
  
  saveRDS(pm_list, paste0('pm_list_ic_', major_weight, '.rds'))
}

#index_list used for fast subsetting
mat <- sparseMatrix(
  i = match(mh$pmid, unique(mh$pmid)),
  j = match(mh$muid, unique(mh$muid)),
  x = 1,
  dimnames = list(unique(mh$pmid), unique(mh$muid))
)
index_list <- as.matrix(mat!=0) %>% split(., seq_len(nrow(.)))
n <- length(index_list)
pmids <- rownames(mat) %>% as.integer()

#Clean up and prepare the loop
rm(mh, mat)

#Soft Cosine similarity function
fnSc <- function(i) {
  focal <- pm_list[[i]]
  focal_mask <- index_list[[i]]
  
  # Boolean mask for where either focal or alter has non-zero entry. Used for
  # subsetting vectors and term-similarity matrix, so we can multiply them in
  # smaller dimension. 
  masks <- lapply(index_list[(i+1):n], \(x) {focal_mask | x})
  
  sc <- numeric(n-i)
  
  for (j in (i+1):n) {
    
    mask <- masks[[j-i]]
    
    #Subset vectors and similarity matrix
    a <- focal[mask]
    b <- pm_list[[j]][mask]
    S <- dm[mask, mask]
    
    #Calculate soft cosine similarity
    sc[j-i] <- (a %*% S %*% b) / (sqrt(a %*% S %*% a) * sqrt(b %*% S %*% b))
  }
  
  return(sc)
}

plan(multisession, workers = 6) 
options(future.globals.maxSize = 4000 * 1024^2)
matrices <- c(
  paste0('dm_ic_sim', 1:5),
  paste0('dm_1_sim', 1:5)
)

# Compute soft cosine depending on 1) term-similarity matrix used, 2) weight
# given to major terms, 3) whether weighted by information content or not.
for (matrix in matrices) {
  
  dm <- readRDS(paste0(matrix, '.rds'))
  
  prefix <- str_remove(matrix, 'dm_')
  
  for (major_weight in 1:10) { fnLoop(pref = matrix, iterator = major_weight)
    
    for (element in c('ic', '1')) {
      
      pm_list <- readRDS(paste0('pm_list_', element, '_', major_weight, '.rds'))
      scs <- future_lapply(1:(n-1), fnSc)
      saveRDS(
        scs,
        paste0('sc_', prefix, '_w', major_weight, '_', element, '.rds')
      )
      rm(scs)
    }
  }
}

plan(sequential)

### 5.3.2 Combine ###############################

#Create vector with filenames
sc_files <- with(
  expand.grid(l = 1:5, e = c('ic', '1'), w = 1:10, dt = c('ic', '1')), 
  paste0('sc_', e, '_sim', l, '_w', w, '_', dt, '.rds')
)

#Convert lists to data.tables
for (sc in sc_files){ fnLoop(sc)
  
  file <- readRDS(sc)
  if (any(sapply(file, function(x) any(is.na(x))))) stop('NA')
  
  df <- rbindlist(lapply(seq_along(file), \(i) {
    data.table(
      pmid_A = pmids[i],
      pmid_B = pmids[(i + 1):length(pmids)],
      sim = file[[i]]
    )
  }))
  
  setnames(df, 'sim', str_remove_all(sc, 'sc_|\\.rds'))
  setkey(df, 'pmid_A', 'pmid_B')
  saveRDS(df, paste0('df_', sc))
  rm(file, df)
}

#Merge data.tables. Done in two rounds due to file sizes
sc_dfs <- paste0('df_', sc_files) %>% 
  split(cut(seq_along(.), 2))

df <- readRDS(sc_dfs[[1]][1])
for (i in sc_dfs[[1]][-1]) { fnLoop(i)
  df <- merge(df, readRDS(i), by = c('pmid_A', 'pmid_B'))
}
saveRDS(df, 'sc_part1.rds')

df <- readRDS(sc_dfs[[2]][1])
for (i in sc_dfs[[2]][-1]) { fnLoop(i)
  df <- merge(df, readRDS(i), by = c('pmid_A', 'pmid_B'))
}
saveRDS(df, 'sc_part2.rds')
rm(df)

df <- merge(
  readRDS('sc_part1.rds'),
  readRDS('sc_part2.rds'),
  by = c('pmid_A', 'pmid_B')
)

saveRDS(df, 'sc.rds')

# Add topics and relevance judgement
load('topics.rda')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>% filter(pmid %in% pmidsWithMh)

bm_data <- lapply(seq_along(topics), fnRjs) %>% bind_rows()
rm(df)
saveRDS(bm_data, file = 'sc_bm.rds')

sc_02 <- bm_data %>% filter(rj1 == 0 & rj2 == 2)
saveRDS(sc_02, 'sc_02.rds')

sc_22 <- bm_data %>% filter(rj1 == 2 & rj2 == 2)
saveRDS(sc_22, 'sc_22.rds')

rm(bm_data, sc_02, sc_22, dm, rjs, sc_files, pmidsWithMh, topics, fnRjs, index_list, pm_list, sc_dfs, element, i, major_weight, matrices, matrix, n, pmids, prefix, sc, t1, t2, fnSc)

# **************************************************************************** #
#******************************************************************************#
# 6. Accuracy test                                                          ####
#******************************************************************************#
# **************************************************************************** #

# ********************************************* #
## 6.1 Test 1 + Summary stats                ####
# ********************************************* #

plan(multisession, workers = 6) 
options(future.globals.maxSize = 2000 * 1024^2)

#Conventional cosine

for (type in c('cc', 'ts', 'sc')) { fnLoop(type)
  
  df02 <- readRDS(paste0(type, '_02.rds')) %>% data.frame()
  df22 <- readRDS(paste0(type, '_22.rds')) %>% data.frame()
  
  vars <- setdiff(colnames(df02), c('pmid_A', 'pmid_B', 'topic', 'rj1', 'rj2'))
  
  clfd <- future_lapply(
    vars,
    \(x) cliff.delta(df22[[x]], df02[[x]])$estimate
  ) %>% unlist()
  
  test1 <- data.frame(
    var      = vars,
    mean02   = colMeans(df02[vars]),
    mean22   = colMeans(df22[vars]),
    median02 = apply(df02[vars], 2, median),
    median22 = apply(df22[vars], 2, median),
    skew02   = apply(df02[vars], 2, skewness, type = 1),
    skew22   = apply(df22[vars], 2, skewness, type = 1),
    sd02     = apply(df02[vars], 2, sd),
    sd22     = apply(df22[vars], 2, sd),
    clfd
  ) %>% 
    mutate(var = str_remove(var, '^X'))
  
  saveRDS(test1, paste0('test1_', type, '.rds'))
}

rm(df02, df22, clfd, test1, type, vars)
plan(sequential)

# ********************************************* #
## 6.3 Test 2                               #####
# ********************************************* #

load('topics.rda')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>%
  filter(pmid %in% pmidsWithMh & rj %in% c(0, 2))

for (type in c('cc', 'sc', 'ts')) { fnLoop(type) 
  
  bm_data <- readRDS(paste0(type, '_bm.rds'))
  
  vars <- setdiff(
    colnames(bm_data), 
    c('pmid_A', 'pmid_B', 'topic', 'rj1', 'rj2')
  )
  
  results_all <- list()
  
  for (j in seq_along(topics)) { fnLoop(iterator = 'topic', suf = j)
    
    topic_number <- topics[j]
    
    df <- bm_data %>% filter(topic == topic_number)
    
    pmids <- union(df$pmid_A, df$pmid_B) 
    
    rjs_topic <- rjs %>% filter(topic == topic_number & pmid %in% pmids)
    topic_pmids <- rjs_topic %>% filter(rj==2) %>% pull(pmid)
    not_topic_pmids <- rjs_topic %>% filter(rj==0) %>% pull(pmid)
    
    results_topic <- list()
    
    for (i in 1:50) { cat(i, ' ')
      
      set.seed(i)
      topic_sample <- topic_pmids %>% sample(size = 10)
      
      set.seed(i)
      not_topic_sample <- not_topic_pmids %>% sample(size = 10)
      
      consider_pmids <- rjs_topic$pmid[
        !rjs_topic$pmid %in% c(topic_sample, not_topic_sample)
      ] 
      
      results_topic[[i]] <- lapply(consider_pmids, function(consider) {
        
        test <- subset(df, pmid_A == consider | pmid_B == consider)
        topic <- subset(
          test, 
          pmid_A %in% topic_sample | pmid_B %in% topic_sample
        )
        not_topic <- subset(
          test, 
          pmid_A %in% not_topic_sample | pmid_B %in% not_topic_sample
        ) 
        
        result <- cbind(
          var = vars,
          truth = rjs_topic$rj[rjs_topic$pmid == consider],
          decision = sapply(vars, function(var) {
            ifelse(max(not_topic[[var]]) < max(topic[[var]]), 2, 0)
          })
        )
        
        return(result)
        
      })
    }
    
    results_all[[j]] <- lapply(
      results_topic,
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
  
  test2 <- results_all %>% bind_rows() %>% 
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
  
  saveRDS(test2, paste0('test2_', type, '.rds'))
  
}

rm(topic_number, df, pmids, rjs_topic, topic_pmids, not_topic_pmids, topic_sample, not_topic_sample, consider_pmids, results_topic, results_all, test2, i, j, rjs, topics, type, vars, bm_data, pmidsWithMh)

# **************************************************************************** #
#******************************************************************************#
# 7. Results                                                                ####
#******************************************************************************#
# **************************************************************************** #

# ********************************************* #
## 7.1 Parameters                            ####
# ********************************************* #

#Function to create plotly heatplot, was used in the explorative phase.
fnHeatPlot <- function(data, performance, tit, xlabs, ylabs = 1:10) {
  
  pdata <- data %>% select(contains(performance)) %>% as.matrix()
  
  l <- length(xlabs)
  
  p <- plot_ly(
    z = pdata, type = 'heatmap', y = 1:10, showscale = F, 
    colorscale = list(c(0, '#fca636'), c(0.25, '#e16462'), c(0.5, '#b12a90'), 
                      c(0.75, '#6a00a8'), c(1, '#0d0887')),
    text = pdata, texttemplate = "%{text}"
  ) %>% 
    layout(
      title = list(text = tit, yanchor = 'top', y = '0.99'),
      yaxis = list(title = 'Major weight', ticktext = ylabs, tickvals = 1:length(ylabs)),
      xaxis = list(tickvals = 0:l, ticktext = xlabs)
    )
  
  return(p)
}

#Common settings for plots in paper
set_com <- list(
  geom_tile(),
  scale_fill_viridis(discrete = FALSE, direction = -1, option = 'A'),
  scale_color_manual(values = c("white", "grey40")),
  theme_minimal(),
  theme(
    legend.position = 'none',
    panel.grid = element_blank(),
    axis.text.y = element_text(margin = margin(r = -7)),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(margin = margin(t = -8)),
    plot.title = element_text(size = 12, hjust = 0.5, margin = margin(b = -6))
  )
)

### 7.1.1 Conventional Cosine ###################

#Data
cc_long <- inner_join(
  readRDS('test1_cc.rds'), readRDS('test2_cc.rds'), by = 'var'
) %>% 
  mutate(
    name = str_extract(var, "[^_]+_[^_]+$"),
    across(where(is.numeric), \(x) round(x, 5))
  ) %>%
  rename_with(.fn = tolower) %>% 
  separate(name, into = c('coord', 'w'), sep = '_') %>% 
  mutate(w = str_extract(w, '\\d+') %>% as.integer()) %>% 
  arrange(coord)
saveRDS(cc_long, 'cc_long.rds')

# Figure in paper
df <- cc_long %>% select(w, var, clfd, phi) %>% 
  mutate(
    vers = str_extract(var, '_1_|_ic_w|_q_'),
    across(where(is.numeric), ~ round(.x, 4))
  )

set_cc <- list(
  set_com,
  scale_x_discrete(labels = c('Basic', 'IC weighted', 'IC weighted\nand qualifiers')),
  scale_y_continuous(breaks = 1:10),
  theme(axis.text.x = element_text(
    angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = -6)
  )),
  labs(x = NULL)
)

cc_win_clfd <- ggplot(df, aes(vers, w, fill = clfd)) +
  set_cc +
  geom_text(aes(
    label = as.character(clfd) %>% str_replace('0.', '.'), color = clfd < 0.35),
    size = 3
  ) +
  labs(y = 'Major weight', title = "Cliff's δ")

cc_win_phi <- ggplot(df, aes(vers, w, fill = phi)) +
  set_cc +
  geom_text(aes(
    label = as.character(phi) %>% str_replace('0.', '.'), color = phi < 0.19), 
    size = 3
  ) +
  labs(y = NULL, title = 'φ')

cc_win_comb <- cc_win_clfd + cc_win_phi
cc_win_comb

ggsave("cc_win.tiff", plot = cc_win_comb, device = "tiff", dpi = 1200, 
       width = 174, height = 80, units = "mm", compression = "lzw")

#Clean
rm(cc, cc_clfd, cc_phi, p, cc_long, cc_win_clfd, cc_win_comb, cc_win_phi, df, set_cc)

### 7.1.2 Soft Cosine ###########################

#Data
sc_long <- inner_join(
  readRDS('test1_sc.rds'), readRDS('test2_sc.rds'), by = 'var'
) %>% 
  mutate(var = str_remove(var, 'f_')) %>% 
  rename_with(.fn = tolower) %>% 
  separate(var, into = c('dt', 'lambda', 'w', 'c'), sep = '_', remove = F) %>%
  filter(!is.na(c)) %>% 
  mutate(
    across(c(phi, clfd), ~ round(.x, 5)),
    lambda = str_extract(lambda, '\\d') %>% as.integer(),
    w = str_extract(w, '\\d+') %>% as.integer()
  )
saveRDS(sc_long, 'sc_long.rds')

#Plot in paper
df <- sc_long %>% filter(dt == 'ic' & c == 'ic') %>% 
  mutate(across(where(is.numeric), ~ round(.x, 4)))

set_sc <- list(
  set_com,
  scale_y_continuous(breaks = 1:10),
  scale_x_continuous(breaks = 1:5),
  labs(x = 'λ')
)

sc_win_clfd <- ggplot(df, aes(lambda, w, fill = clfd)) +
  set_sc +
  geom_text(aes(
    label = as.character(clfd) %>% str_replace('0.', '.'), color = clfd < 0.28), 
    size = 3
  ) +
  labs(title = "Cliff's δ", y = 'Major weight')

sc_win_phi <- ggplot(df, aes(lambda, w, fill = phi)) +
  set_sc +
  geom_text(aes(
    label = as.character(phi) %>% str_replace('0.', '.'), color = phi < 0.165), 
    size = 3
  ) +
  labs(title = "φ", y = NULL)

sc_win_comb <- sc_win_clfd + sc_win_phi +
  plot_layout(guides = "collect", axes = "collect_x")
sc_win_comb

ggsave("sc_win.tiff", plot = sc_win_comb, device = "tiff", dpi = 1200, 
       width = 174, height = 75, units = "mm", compression = "lzw")

#Plots which were used in exploratory phase
sc <- sc_long %>% select(dt, lambda, w, c, clfd, phi)

fnSc <- function(dt_val, c_val) {
  t1 <- sc %>% 
    filter(dt == dt_val & c == c_val) %>% 
    select(lambda, w, clfd, phi) %>% 
    pivot_wider(names_from = lambda, values_from = c(clfd, phi))
}


dt1_c1 <- fnSc('1', '1')
dt1_cic <- fnSc('1', 'ic')
dtic_c1 <- fnSc('ic', '1')
dtic_cic <- fnSc('ic', 'ic')

# Since the plots have a lot of common settings, we create a function to capture
# these. Same is done for maximum term similarities. 
fnScPlot <- function(df, pfmc, dt, c) {
  p <- fnHeatPlot(
    data = df, performance = pfmc,
    tit = paste0(
      "Cliff's d, Soft cosine, Distance-terms = ", dt, 
      ' & Coordinates = ', c
    ),
    xlabs = 1:5
  )
}

# Naming:
#sc: These plots concern soft cosine. 
#dt: Dist terms (how is it measured: delta-IC or just 1)
#c : coordinates (of the vectors, are they IC or 1)
#last underscore: is the performance measure clfd (Cliff's D) or phi 

sc_dt1_c1_clfd <- fnScPlot(dt1_c1, 'clfd', '1', '1')
sc_dt1_cic_clfd <- fnScPlot(dt1_cic, 'clfd', '1', 'IC')
sc_dtic_c1_clfd <- fnScPlot(dtic_c1, 'clfd', 'IC', '1')
sc_dtic_cic_clfd <- fnScPlot(dtic_cic, 'clfd', 'IC', 'IC')

sc_dt1_c1_phi <- fnScPlot(dt1_c1, 'phi', '1', '1')
sc_dt1_cic_phi <- fnScPlot(dt1_cic, 'phi', '1', 'IC')
sc_dtic_c1_phi <- fnScPlot(dtic_c1, 'phi', 'IC', '1')
sc_dtic_cic_phi <- fnScPlot(dtic_cic, 'phi', 'IC', 'IC')

for (p in c('sc_dtic_c1_clfd', 'sc_dtic_cic_clfd', 'sc_dt1_c1_clfd', 'sc_dt1_cic_clfd', 'sc_dtic_c1_phi', 'sc_dtic_cic_phi', 'sc_dt1_c1_phi', 'sc_dt1_cic_phi')) {
  saveWidget(get(p), file = paste0(p, '.html'), selfcontained = FALSE)
}

#Clean
rm(sc_dtic_c1_clfd, sc_dtic_cic_clfd, sc_dt1_c1_clfd, sc_dt1_cic_clfd, sc_dtic_c1_phi, sc_dtic_cic_phi, sc_dt1_c1_phi, sc_dt1_cic_phi, sc, fnSc, p, dt1_c1, dt1_cic, dtic_c1, dtic_cic, fnScPlot, sc_long, df, sc_win_clfd, sc_win_comb, sc_win_phi, set_sc)

### 7.1.3 Term similarity #######################

ts_long <- inner_join(
  readRDS('test1_ts.rds'), readRDS('test2_ts.rds'), by = 'var'
) %>% 
  rename_with(.fn = tolower) %>% 
  mutate(
    across(c(phi, clfd), ~ round(.x, 5)),
    var = ifelse(str_detect(var, 's_'), paste0(var, '_w0'), var),
    w = str_extract(var, '(?<=w)\\d+') %>% as.integer(),
    lambda = str_extract(var, '(?<=_(sim))\\d') %>% as.integer(),
    lambda = ifelse(is.na(lambda), 0, lambda),
    dt = case_when(
      str_detect(var, '_1_') ~ '1',
      str_detect(var, '_ic_') ~ 'ic'
    )
  ) %>% 
  arrange(w)
saveRDS(ts_long, 'ts_long.rds')

# Figure in paper
df <- ts_long %>% filter(dt == 'ic' & str_detect(var, '_ic_')) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 4)))

set_ts <- list(
  set_com,
  scale_y_continuous(breaks = 0:20, labels = c("Slim", 1:20)),
  scale_x_continuous(breaks = 0:5, labels = c('Distance', paste0('λ', 1:5))),
  labs(x = NULL)
)

ts_win_clfd <- ggplot(df, aes(lambda, w, fill = clfd)) +
  set_ts +
  geom_text(aes(
    label = as.character(clfd) %>% str_replace('0.', '.'), color = clfd < 0.29),
    size = 3
  ) +
  labs(y = 'Major weight', title = "Cliff's δ")
ts_win_clfd

ts_win_phi <- ggplot(df, aes(lambda, w, fill = phi)) +
  set_ts +
  geom_text(aes(
    label = as.character(phi) %>% str_replace('0.', '.'), color = phi < 0.15),
    size = 3
  ) +
  labs(y = NULL, title = "φ")
ts_win_phi

ts_win_comb <- ts_win_clfd + ts_win_phi
ts_win_comb

ggsave("ts_win.tiff", plot = ts_win_comb, device = "tiff", dpi = 1200, 
       width = 174, height = 110, units = "mm", compression = "lzw")

#Plots used in explorative phase
ts <- ts_long %>% select(var, dt, w, lambda, clfd, phi)

fnTs <- function(dt_val) {
  ts %>% filter(dt == dt_val) %>% 
    select(w, lambda, clfd, phi) %>% 
    pivot_wider(names_from = lambda, values_from = c(clfd, phi))
}

ts_ic <- fnTs('ic')
ts_1 <- fnTs('1')
ts_wp <- fnTs('wp')
ts_jc <- fnTs('jc')

fnTsPlot <- function(df, pfmc, dt) {
  
  p <- fnHeatPlot(
    data = df, performance = pfmc,
    xlabs = c('Distance', 1:5), ylabs = c('slim', 1:10), 
    tit = paste0(
      ifelse(pfmc == 'clfd', "Cliff's d", 'Phi'), 
      ', Term similarity, Distance-terms = ',
      dt
    )
  ) %>% layout(xaxis = list(title = 'Lambda'))
}

#Plots. Cliff's d and phi, then to .html. 
ts_dtic_clfd  <- fnTsPlot(ts_ic, 'clfd', 'ic')
ts_dt1_clfd   <- fnTsPlot(ts_1, 'clfd', '1')

ts_dtic_phi  <- fnTsPlot(ts_ic, 'phi', 'ic')
ts_dt1_phi   <- fnTsPlot(ts_1, 'phi', '1')

for (p in c('ts_dtic_phi', 'ts_dt1_phi', 'ts_dtjc_phi', 'ts_dtic_clfd', 'ts_dt1_clfd')) {
  saveWidget(get(p), paste0(p, '.html'), selfcontained = F)
}

#Clean
rm(ts, ts_ic, ts_1, ts_jc, ts_wp, ts_dtic_phi, ts_dt1_phi, ts_dtic_clfd, ts_dt1_clfd, p, fnTs, fnTsPlot, ts_long, fnHeatPlot, ts_win_clfd, ts_win_phi, ts_win_comb, set_ts, df)

# ********************************************* #
## 7.2 Table 1                               ####
# ********************************************* #

cc <- readRDS('cc_long.rds') %>% mutate(comb = phi + clfd)
sc <- readRDS('sc_long.rds') %>% mutate(comb = phi + clfd)
ts <- readRDS('ts_long.rds') %>% mutate(comb = phi + clfd)
conv <- readRDS('cc_long.rds') %>% filter(var == 'cc_1_w1') %>% 
  select(mean02, mean22, skew02, skew22, var, w, clfd, phi)

#Winners
cc_win <- cc %>% group_by(coord) %>% filter(comb == max(comb)) %>% 
  rename(c = coord) 
sc_win <- sc %>% group_by(dt, c) %>% filter(comb == max(comb)) 
ts_win <- ts %>% group_by(dt) %>% filter(comb == max(comb))

#Table 1
tab1 <- bind_rows(conv, cc_win, sc_win, ts_win) %>% 
  select(var, c, dt, w, lambda, clfd, phi, mean02, skew02, mean22, skew22) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))

rm(cc, cc_win, conv, sc, sc_win, tab1, ts, ts_win)

# ********************************************* #
## 7.3 Distribution                          ####
# ********************************************* #

### Create data **********************#

fnDist <- function(cols, type) {
  colsFn <- c('pmid_A', 'pmid_B', cols)
  p02 <- paste0(type, '_02.rds')
  p22 <- paste0(type, '_22.rds')
  
  df <- bind_rows(
    readRDS(p02) %>% select(all_of(colsFn)) %>% mutate(same = F),
    readRDS(p22) %>% select(all_of(colsFn)) %>% mutate(same = T),
  ) %>% 
    mutate(
      pmid_a = pmin(pmid_A, pmid_B),
      pmid_b = pmax(pmid_A, pmid_B)
    ) %>% 
    select(-pmid_A, -pmid_B) %>% 
    arrange(pmid_a, pmid_b)
}

cc <- fnDist(c('cc_ic_w2', 'cc_1_w1'), 'cc')
sc <- fnDist(c('ic_sim1_w3_ic'),       'sc')
ts <- fnDist(c('f_ic_sim2_w16'),       'ts')

keys <- c('pmid_a', 'pmid_b', 'same')
dist <- merge(cc, sc, by = keys) %>% merge(ts, by = keys) %>% 
  select(-pmid_a, -pmid_b) %>% 
  pivot_longer(cols = -same, names_to = 'var', values_to = 'val')

rm(cc, sc, ts, keys, fnDist)
saveRDS(dist, 'density_distribution.rds')

### Common settings ******************#

dist <- readRDS('density_distribution.rds')

order <- c('ic_sim1_w3_ic', 'cc_ic_w2', 'f_ic_sim2_w16', 'cc_1_w1')
labs <- c(
  'ic_sim1_w3_ic' = "IC soft cosine w3 λ1",
  'cc_ic_w2' = "IC Salton's Cosine w2",
  'f_ic_sim2_w16' = "IC maximum term similarities w16 λ2",
  'cc_1_w1' = "Standard Salton's cosine"
)
lines <- c(
  'ic_sim1_w3_ic' = "solid",
  'cc_ic_w2' = "twodash",
  'f_ic_sim2_w16' = "longdash",
  'cc_1_w1' = "dotted"
)

set_com <- list(
  geom_density(),
  scale_color_viridis_d(
    name = "Method", 
    labels = labs, 
    limits = order, 
    option = "D", 
    direction = -1
  ),
  scale_linetype_manual(
    name = "Method",
    values = lines,
    labels = labs,
    limits = order
  ),
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 0.25, 0.5, 0.75)),
  coord_cartesian(xlim = c(0, 0.75)),
  theme_minimal(),
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = 'grey70'),
    axis.ticks = element_line(color = 'grey70'),
    axis.ticks.length = unit(0.15, "cm"),
    plot.title = element_text(size = 12)
  ),
  labs(x = 'Relatedness', y = NULL)
)

### Create plots *********************#

#Left plot (seperate)
p_sep <- ggplot(dist[!dist$same, ], aes(x = val, color = var, linetype = var)) +
  set_com +
  scale_y_continuous(
    breaks = seq(0, 18, by = 2), limits = c(0, 18), expand = c(0, 0)
  ) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
  ) +
  labs(color = 'Method', linetype = 'Method', title = 'Seperate topics')
p_sep

#Right plot (same)
p_sam <- ggplot(dist[dist$same, ], aes(x = val, color = var, linetype = var)) +
  set_com +
  scale_y_continuous(breaks = seq(0, 6, by = 2), expand = c(0, 0)) +
  theme(legend.position = 'none') +
  labs(title = 'Same topics')

p_sam

comb <- p_sep + p_sam + plot_layout(axis_titles = 'collect')
comb
ggsave("distr_samXsep.tiff", plot = comb, device = "tiff", dpi = 1200, 
       width = 174, height = 90, units = "mm", compression = "lzw")

rm(comb, dist, p_sam, p_sep, set_com, labs, lines, order)

# **************************************************************************** #
#******************************************************************************#
# 8. Misc stats                                                             ####
#******************************************************************************#
# **************************************************************************** #

#Number of publications + % relevant judgments
load('topics.rda')
load('pmidsWithMh.rda')
rjs <- read.csv('rjs6.txt') %>%
  filter(pmid %in% pmidsWithMh & rj %in% c(0, 2))

used <- rjs %>% filter(topic %in% topics)

length(used$pmid %>% unique())

(sum(used$rj == 2) / nrow(used)) %>% round(3)

#Cliff's d: Number of publication pairs + % relevant-relevant pairs
nr_r <- readRDS('cc_02.rds') %>% nrow()
r_r <- readRDS('cc_22.rds') %>% nrow()

n <- nr_r + r_r
(r_r/n) %>% round(3)

#Phi: Number of decisions.
n <- readRDS('test2_cc.rds')
n$n

rm(used, nr_r, r_r, n, rjs, pmidsWithMh, topics)
