rm(list=ls(all=T))

# set options
options(stringsAsFactors = F)


#################
# PRELIMINARIES #
#################

# install and load packages
sapply(c("dplyr", "devtools", "koRpus", "data.table", "reshape2", "pbapply"), 
       function(x) if(!is.element(x, installed.packages())) install.packages(x, dependencies = T))


# collostructions
if(!is.element("collostructions", installed.packages())) {
  if(.Platform$OS.type!="Windows") {
    install.packages("http://userpage.fu-berlin.de/~flach/wp-content/uploads/collostructions_0.0.10.tar.gz", repos = NULL)
  } else {
    install.packages("http://userpage.fu-berlin.de/~flach/wp-content/uploads/collostructions_0.0.10.zip", repos = NULL)
  }
}

# concordances
if(!is.element("concordances", installed.packages())) {
  devtools::install_github("hartmast/concordances")
}

lapply(list("dplyr", "collostructions", "concordances", "koRpus", "data.table", "pbapply"), 
       require, character.only=T)


# read in files
can <- getNSE("ngrams_encow/can.txt")
gb <- getNSE("ngrams_encow/gb.txt")
us <- getNSE("ngrams_encow/us.txt")


# make n-gram frequency lists
can1 <- can$Tag2 %>% table %>% as.data.frame()
colnames(can1) <- c("Lemma", "Freq_can")

gb1 <- gb$Tag2 %>% table %>% as.data.frame()
colnames(gb1) <- c("Lemma", "Freq_gb")

us1 <- us$Tag2 %>% table %>% as.data.frame()
colnames(us1) <- c("Lemma", "Freq_us")

#combine frequency lists
fl <- merge(merge(gb1, can1, all = T), us1, all=T)
fl[is.na(fl)] <- 0


# collexeme analysis
collex.dist(fl, threshold = 10) # not very informative...


############################
# Toy example: unmitigated #
############################

# read concordance "unmitigated" + x
unm <- getCWB("unmitigated/unmitigated.txt")

# read BNC lemma frequency list
lemmas_bnc <- read.table("lemmalist_bnc.txt", fill=T, quote="")
colnames(lemmas_bnc) <- c("Freq", "Lemma", "whatever")

# get collexemes of unmitigated: replace "unmitigated" by nothing...
unm$tagB <- gsub("unmitigated ", "", unm$tagB) #lemmas are in tagB

# make table
unm_t <- unm$tagB %>% table %>% as.data.frame
colnames(unm_t) <- c("Lemma", "Freq")

# add corpus frequencies
unm_t$cf <- NA

for(i in 1:nrow(unm_t)) {
  if(unm_t$Lemma[i] %in% lemmas_bnc$Lemma) {
    unm_t$cf[i] <- sum(lemmas_bnc[which(lemmas_bnc$Lemma==unm_t$Lemma[i]),]$Freq)
  } else {
    unm_t$cf[i] <- 0
  }
}


# collex - corpsize = frequency of nouns, assessed
# via CQP query [pos="N.*"] --> count Last by pos

collex(unm_t, corpsize = 14430588 + 5218330)



###############################
# Toy example: into-causative #
###############################

into_coca <- read.table("into-causative/into_coca.csv", sep="\t", head=T, 
                        quote="", fill=T, stringsAsFactors = F)
into_bnc <- getCWB("into-causative/into_bnc.txt")
head(into_bnc) # tokens are in column tagA, lemmas (hw) in columns tagB


## AMERICAN ENGLISH DATA

# The BYU data have to be prepared a bit:
# we need the *last* word of the left context (= the verb)
# and the last word of the keyword column.

into_coca$progressive <- into_coca$verb <- NA

for(i in 1:nrow(into_coca)) {
  x <- unlist(strsplit(into_coca$left[i], " "))[length(unlist(strsplit(into_coca$left[i], " ")))]
  if(length(x)>0) {
    into_coca$verb[i] <- tolower(gsub("[[:space:]]", "", x))
  }
  
  y <- unlist(strsplit(into_coca$key[i], " "))[length(unlist(strsplit(into_coca$key[i], " ")))]
  
  if(length(y)>0){
    into_coca$progressive[i] <- tolower(gsub("[[:space:]]", "", y))
  }
  
}


# some columns are empty, we have to remove them:
into_coca <- into_coca[which(gsub("[[:punct:]]", "", into_coca$verb)!=""),]
into_coca <- into_coca[which(gsub("[[:punct:]]", "", into_coca$progressive)!=""),]


# unfortunately, we can't export lemmas from the BYU interface.
# So we just lemmatize the vector ourselves.
# Be sure to change the path to your installation of the TreeTagger:
into_coca$verb_lemma <- treetag(into_coca$verb, treetagger="/Users/stefanhartmann/Downloads/TreeTagger/cmd/tree-tagger-english",
        lang = "en", format = "obj")@TT.res$lemma
 
# covarying collexeme analysis
into_coca %>% select(verb_lemma, progressive) %>% collex.covar()



## BRITISH ENGLISH DATA
into_bnc$verb_lemma <- into_bnc$progressive <- NA

into_bnc$verb_lemma <- sapply(1:nrow(into_bnc), function(i)
  unlist(strsplit(unlist(strsplit(into_bnc$Keyword[i], " "))[1], "/"))[2])
into_bnc$progressive <- sapply(1:nrow(into_bnc), function(i)
  unlist(strsplit(unlist(strsplit(into_bnc$Keyword[i], " "))[length(unlist(strsplit(into_bnc$Keyword[i], " ")))], "/"))[1])


# distinctive collexeme analysis
into_bnc %>% select(verb_lemma, progressive) %>% collex.covar()


# verbs: British vs. American
brit <- into_bnc$verb_lemma %>% table %>% as.data.frame
colnames(brit) <- c("Lemma", "Freq_brit")

am <- into_coca$verb_lemma %>% table %>% as.data.frame()
colnames(am) <- c("Lemma", "Freq_am")

# combine dfs:
both <- merge(brit, am, all = T)
both[is.na(both)] <- 0
both %>% collex.dist(reverse = F)


################################
# Toy example: n-gram analyses #
################################


## UNIGRAM ANALYSIS MALE VS. FEMALE ##
female <- fread("https://onedrive.live.com/download?cid=D5FA76D941B63BA2&resid=D5FA76D941B63BA2%2146750&authkey=ALDO-hdF-qkZrQU", header = F)
male <- fread("https://onedrive.live.com/download?cid=D5FA76D941B63BA2&resid=D5FA76D941B63BA2%2146751&authkey=AAzlmECWFSgJNLg", header = F)

female$V2 <- pbsapply(1:nrow(female), function(i) gsub(">", "", unlist(strsplit(as.character(female$V2)[i], "/"))[3]))
male$V2 <- pbsapply(1:nrow(male), function(i) gsub(">", "", unlist(strsplit(as.character(male$V2)[i], "/"))[3]))


# make tables
fem <- female[,2] %>% table %>% as.data.frame()
colnames(fem) <- c("word", "freq_fem")
ml <- male[,2] %>% table %>% as.data.frame()
colnames(ml) <- c("word", "freq_m")

# merge tables
both <- merge(fem, ml, all = T)

# replace NAs by 0
both[is.na(both)] <- 0

both %>% collex.dist(reverse = T, am = "chisq")


