
#------------------------------------------------------------
# similarity measure
#------------------------------------------------------------

similarity <- function(seqA, seqB) {
  # levenshtein distance between seqA and seqB
  seqdistance <- stringdist(seqA, seqB, "lv")
  # similarity score between seqA and seqB
  # normalise levdist using the length of the longest sequence
  longest <- tibble(lenA = nchar(seqA), lenB = nchar(seqB)) %>% 
    rowwise() %>% 
    mutate(longest = max(lenA, lenB)) %>% 
    pull(longest)
  # longest <- max(nchar(seqA), nchar(seqB))
  seqsimilarity <- 1 - (seqdistance / longest)
  return(seqsimilarity)
}

#------------------------------------------------------------
# summarySE
#------------------------------------------------------------

# partially based on: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence 
## interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)

summarySE <- function(data, measurevar, groupvars, conf.interval=.95) {
  
  # Mean by groups
  datac <- data %>% 
    mutate_(measurevar = measurevar) %>% 
    group_by_(.dots = groupvars) %>%
    summarise(mean = mean(measurevar, na.rm = T), sd = sd(measurevar, na.rm = T), N = n()) %>% 
    ungroup()
  datac[,measurevar] <- datac$mean
  
  # Calculate standard error of the mean
  datac$se <- datac$sd / sqrt(datac$N)
  
  # Confidence interval multiplier for standard error
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#--------------------------------------------------------------
# combis()
#--------------------------------------------------------------

# all possible grams of size i (e.g. 1grams to 5grams)
combis <- function(min.n, max.n, symbs, clps = "") {
  pgrams <- tibble()
  for (i in min.n:max.n) {
    pgrams1 <- expand.grid(rep(list(symbs), i), stringsAsFactors = F)
    pgrams2 <- apply(pgrams1, 1, paste, collapse = clps)
    pgrams3 <- tibble(gramsize = i, ngram = pgrams2)
    pgrams <- bind_rows(pgrams, pgrams3)
  }
  return(pgrams)
}
