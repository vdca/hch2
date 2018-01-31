
#--------------------------------------------------------
# Globals
#--------------------------------------------------------

rm(list=ls())     # remove previous objects from workspace

source("hch2_funs.R")

library(tidyverse)
theme_set(theme_bw() + theme(strip.background = element_blank()))
library(broom)
library(lme4)
library(lmerTest)
library(stringdist)

#--------------------------------------------------------
# Read data
#--------------------------------------------------------

# all sequences in all 11 (initial random+10) generations
d <- readRDS("hch2_sequences.rds") %>% tbl_df() %>%
  select(chain, generation, subjID, seqid = seqid1, seqIN = inseq, seqOUT = outseq)
d %>% write_csv("hch2_sequences.csv")
d <- read_csv("hch2_sequences.csv", col_types = "iicicc") 

#--------------------------------------------------------
#' # Learnability (plots)
#--------------------------------------------------------

#' ** strangely enough, models give problems if score is *un*-rounded
d <- d %>% 
  mutate(score = similarity(seqIN, seqOUT)) %>% 
  mutate(score = round(score, 2))

dh <- d %>% filter(generation != 0)

#' summarise learnability by chain and generation, with confidence intervals,
#' then average across chains
sylscore <- summarySE(dh, measurevar="score",
                      groupvars=c("chain", "generation"))
sylscore2 <- summarySE(sylscore, measurevar="mean", groupvars=c("generation"))

similarity1 <- ggplot(sylscore2) +
  aes(x = as.factor(generation), y = mean, group = 1) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_line() + # geom_point() +
  labs(x = "Generation", y = "Normalised similarity")
similarity1

similarity2 <- sylscore %>% 
  mutate(chain = paste("Chain", chain)) %>% 
  crossing(highlight = unique(.$chain)) %>% 
  ggplot() +
  aes(x = as.factor(generation), y = mean, group = chain, alpha = highlight == chain) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_line() + # geom_point() +
  facet_wrap(~highlight) + 
  scale_alpha_discrete(range = c(.2, 1)) +
  theme(legend.position="none") +
  labs(x = "Generation", y = "Normalised similarity")
similarity2

#--------------------------------------------------------
#' # Learnability (models)
#--------------------------------------------------------

# full
mm.sim.s <- lmer(score ~ generation + (1+generation|chain), dh, REML = F); summary(mm.sim.s)
# null
mm.sim.null.s <- lmer(score ~ 1 + (1+generation|chain), dh, REML = F)
# comparison
anova(mm.sim.null.s, mm.sim.s)

#--------------------------------------------------------
#' # Length (plots)
#--------------------------------------------------------

d <- d %>% 
  mutate(lenIN = nchar(seqIN),
         lenOUT = nchar(seqOUT))
dh <- d %>% filter(generation != 0)

#' summarise sequence length by chain and generation, with confidence intervals,
#' then average across chains

lensum <- summarySE(d, "lenOUT", c("generation", "chain"))
lensum2 <- summarySE(lensum, "lenOUT", c("generation"))

length1 <- ggplot(lensum2) + aes(x = as.factor(generation), y = lenOUT, group = 1) +
  geom_line() +
  geom_errorbar(aes(ymin=lenOUT-ci, ymax=lenOUT+ci), width=.1) +
  labs(x = "Generation", y = "Length of output sequences (n of syllables)")
length1

length2 <- lensum %>% 
  mutate(chain = paste("Chain", chain)) %>% 
  crossing(highlight = unique(.$chain)) %>% 
  ggplot() +
  aes(x = as.factor(generation), y = mean, group = chain, alpha = highlight == chain) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_line() + # geom_point() +
  facet_wrap(~highlight) + 
  scale_alpha_discrete(range = c(.2, 1)) +
  theme(legend.position="none") +
  labs(x = "Generation", y = "Length of output sequences (n of syllables)")
length2

#--------------------------------------------------------
#' # Length (models)
#--------------------------------------------------------

#' Does average length of produced sequences decrease?

# full
mm.len.s <- lmer(lenOUT ~ generation + (1+generation|chain), d, REML = F); summary(mm.len.s)
# null
mm.len.null.s <- lmer(lenOUT ~ 1 + (1+generation|chain), d, REML = F)
# comparison
anova(mm.len.null.s, mm.len.s)

#' Only with human data? (without computer-generated sequences)

# full
mm.len.s.h <- lmer(lenOUT ~ generation + (1+generation|chain), dh, REML = F); summary(mm.len.s.h)
# null
mm.len.null.s.h <- lmer(lenOUT ~ 1 + (1+generation|chain), dh, REML = F)
# comparison
anova(mm.len.null.s.h, mm.len.s.h)

#' Is generation still a significant predictor of learnability
#' after taking into account decreasing length of input sequences?

# full
mm.sim2.s <- lmer(score ~ generation + lenIN + (1+generation+lenIN|chain), dh, REML = F); summary(mm.sim2.s)
# null
mm.sim2.null.s <- lmer(score ~ lenIN + (1+generation+lenIN|chain), dh, REML = F)
# comparison
anova(mm.sim2.null.s, mm.sim2.s)

#--------------------------------------------------------
#' # Set dispersion (plots)
#--------------------------------------------------------

#' compute dispersion within each sequence set.
#' first, generate all pairwise combinations for each sequence set.
disp.combinations <- d %>% 
  tbl_df() %>% 
  select(generation, chain, seqOUT, seqid) %>% 
  mutate(seqOUT = paste(seqid, seqOUT, sep = ".")) %>% 
  group_by(generation, chain) %>% 
  tidyr::expand(seqOUT, seqOUT) %>% 
  ungroup()
#' second, compute dispersion for each sequence pair
disp <- disp.combinations %>% 
  separate(seqOUT, c("seqidA", "outseqA")) %>% 
  separate(seqOUT1, c("seqidB", "outseqB")) %>% 
  mutate(similar = similarity(outseqA, outseqB),
         dispersion = 1 - similar)

#' Dispersion measures how similar each sequence is to other sequences within that set.
#' It's the opposite of similarity. If all the sequences within the set are identical,
#' the overall dispersion score will be 0.

dispersion.chains <- summarySE(disp, "dispersion", c("chain", "generation"))
dispersion <- summarySE(dispersion.chains, "dispersion", "generation")

disp1 <- ggplot(dispersion) + aes(x = as.factor(generation), y = dispersion, group = 1) + 
  geom_errorbar(aes(ymin=dispersion-ci, ymax=dispersion+ci), width=.1) +
  geom_line() +
  labs(x = "Generation", y = "Set dispersion")
disp1

disp2 <- dispersion.chains %>% 
  mutate(chain = paste("Chain", chain)) %>% 
  crossing(highlight = unique(.$chain)) %>% 
  ggplot() +
  aes(x = as.factor(generation), y = mean, group = chain, alpha = highlight == chain) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_line() + # geom_point() +
  facet_wrap(~highlight) + 
  scale_alpha_discrete(range = c(.2, 1)) +
  theme(legend.position="none") +
  labs(x = "Generation", y = "Set dispersion")
disp2

#--------------------------------------------------------
#' # Set dispersion (models)
#--------------------------------------------------------

#' Does set dispersion decrease?

# full
mm.disp <- lmer(dispersion ~ generation + (1|chain), dispersion.chains, REML = F); summary(mm.disp)
# null
mm.disp.null <- lmer(dispersion ~ 1 + (1|chain), dispersion.chains, REML = F)
# comparison
anova(mm.disp.null, mm.disp)

#' Only with human data?

dispersion.chains2 <- dispersion.chains %>% 
  filter(generation != 0)

# full
mm.disp2 <- lmer(dispersion ~ generation + (1|chain), dispersion.chains2, REML = F); summary(mm.disp2)
# null
mm.disp.null2 <- lmer(dispersion ~ 1 + (1|chain), dispersion.chains2, REML = F)
# comparison
anova(mm.disp.null2, mm.disp2)

#--------------------------------------------------------
#' # Dispersion within sequences (plots)
#--------------------------------------------------------

#' dispersion measures increasing similarity (repetition of patterns)
#' within sets of sequences.
#' we now test the increase of repetition within each sequence
#' by measuring the similarity of the first half with the second half
#' of each sequence

d <- d %>% 
  mutate(outseq.1 = str_sub(seqOUT, 0, lenOUT %/% 2),
         outseq.2 = str_sub(seqOUT, (lenOUT %/% 2) + 1, lenOUT),
         sim.1_2 = stringdist(outseq.1, outseq.2, "lv"),
         sim.1_2 = 1-(sim.1_2 / nchar(outseq.2)),
         disp.1_2 = 1 - sim.1_2)
# only human data (without computer-generated sequences)
dh <- d %>% filter(generation != 0)

# carefully control how the split in two halves is done.

# outseq <- "1231234"
# outseq <- "12341234"

# if nchar == odd, then there may be an overlap of 1 symbol.
# this is a problem because the internal dispersion of odd sequences
# would get artificially decreased, and computer sequences are never odd (always 12)
# a1 <- str_sub(outseq, 0, (nchar(outseq)+1) %/% 2); a1
# a2 <- str_sub(outseq, (nchar(outseq)+1) %/% 2, nchar(outseq)); a2

# better to have no overlap between the two halves of the sequence
# and know that the in odd cases the 2nd half will be longer (e.g. 3+4)
# a1 <- str_sub(outseq, 0, nchar(outseq) %/% 2); a1
# a2 <- str_sub(outseq, (nchar(outseq) %/% 2)+1, nchar(outseq)); a2

indisp <- summarySE(d, measurevar="disp.1_2", groupvars=c("chain", "generation"))
indisp2 <- summarySE(indisp, measurevar="mean", groupvars=c("generation"))

indispp1 <- ggplot(indisp2) +
  aes(x = as.factor(generation), y = mean, group = 1) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_line() + # geom_point() +
  labs(x = "Generation", y = "Sequence dispersion")
indispp1

indispp2 <- indisp %>% 
  mutate(chain = paste("Chain", chain)) %>% 
  crossing(highlight = unique(.$chain)) %>% 
  ggplot() +
  aes(x = as.factor(generation), y = mean, group = chain, alpha = highlight == chain) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_line() + # geom_point() +
  facet_wrap(~highlight) + 
  scale_alpha_discrete(range = c(.2, 1)) +
  theme(legend.position="none") +
  labs(x = "Generation", y = "Sequence dispersion")
indispp2

#--------------------------------------------------------
#' # Dispersion within sequences (models)
#--------------------------------------------------------

#' ## Does sequence-internal dispersion decrease?

# full
mm.indisp <- lmer(disp.1_2 ~ generation + (1|chain), d, REML = F); summary(mm.indisp)
# null
mm.indisp.null <- lmer(disp.1_2 ~ 1 + (1|chain), d, REML = F)
# comparison
anova(mm.indisp.null, mm.indisp)

#' Only with human data?

# full
mm.indisp2 <- lmer(disp.1_2 ~ generation + (1|chain), dh, REML = F); summary(mm.indisp2)
# null
mm.indisp.null2 <- lmer(disp.1_2 ~ 1 + (1|chain), dh, REML = F)
# comparison
anova(mm.indisp.null2, mm.indisp2)

#--------------------------------------------------------
#' # Compression (plots)
#--------------------------------------------------------

dcomp <- readRDS("hch2_compression.rds")
dcomp2 <- dcomp %>% filter(generation != 0)

dcsum2 <- summarySE(dcomp, "comp.zlib", "generation")

comp1 <- ggplot(dcsum2) + aes(x = factor(generation), y = mean, group = 1) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  labs(x = "Generation", y = "Compression ratio")
comp1

comp2 <- dcomp %>% 
  mutate(chain = paste("Chain", chain)) %>% 
  crossing(highlight = unique(.$chain)) %>% 
  ggplot() +
  aes(x = as.factor(generation), y = comp.zlib, group = chain, alpha = highlight == chain) + 
  geom_line() + # geom_point() +
  facet_wrap(~highlight) +
  scale_alpha_discrete(range = c(.2, 1)) +
  theme(legend.position="none") +
  labs(x = "Generation", y = "Compression ratio")
comp2

#--------------------------------------------------------
#' # Compression (models)
#--------------------------------------------------------

#' Does compression ratio decrease?

# full
mm.comp <- lmer(comp.zlib ~ generation + (1|chain), dcomp, REML = F); summary(mm.comp)
# null
mm.comp.null <- lmer(comp.zlib ~ 1 + (1|chain), dcomp, REML = F)
# comparison
anova(mm.comp.null, mm.comp)

#' Only with human data?

# full
mm.comp2 <- lmer(comp.zlib ~ generation + (1|chain), dcomp2, REML = F); summary(mm.comp2)
# null
mm.comp.null2 <- lmer(comp.zlib ~ 1 + (1|chain), dcomp2, REML = F)
# comparison
anova(mm.comp.null2, mm.comp2)

#--------------------------------------------------------
#' # Identifiability (plots)
#--------------------------------------------------------

#' ** optimise below loop

# for each sequence, compute similarity with
# same generation, same chain sequences,
# and with same generation, different chain sequences
for (s in unique(d$generation)) {
  for (ch in unique(d$chain)) {
    # all sequences produced by subject
    ds <- d[d$generation == s & d$chain == ch,]
    # all sequences from other chains of same generation
    dother <- d[d$generation == s & d$chain != ch,] %>% pull(seqOUT)
    for (q in unique(ds$seqid)) {
      # all sequences by subject, excluding the reference one
      dwithin <- ds[ds$seqid != q,] %>% pull(seqOUT)
      # reference sequence
      dref <- ds[ds$seqid == q,] %>% pull(seqOUT)
      # compare ref seq to others from same generation/chain
      a <- expand.grid(dref, dwithin, stringsAsFactors = F)
      a$sim.within <- similarity(a$Var1, a$Var2)
      d$sim.within[d$chain == ch & d$generation == s & d$seqid == q] <- mean(a$sim.within)
      # compare ref seq to others from same generation, different chain
      b <- expand.grid(dref, dother, stringsAsFactors = F)
      b$sim.other <- similarity(b$Var1, b$Var2)
      d$sim.other[d$chain == ch & d$generation == s & d$seqid == q] <- mean(b$sim.other)
    }
  }
}

d <- d %>% mutate(ident = sim.within / (sim.within + sim.other))
dh <- d %>% filter(generation != 0)

# compute identifiability score
ident2.d <- summarySE(d, "ident", c("generation", "chain"))
ident1.d <- summarySE(ident2.d, "ident", "generation")

ident1 <- ggplot(ident1.d) + aes(x = as.factor(generation), y = ident, group = 1) + 
  geom_errorbar(aes(ymin=ident-ci, ymax=ident+ci), width=.1) +
  geom_line() +
  labs(x = "Generation", y = "Identifiability")
ident1

ident2 <- ident2.d %>% 
  mutate(chain = paste("Chain", chain)) %>% 
  crossing(highlight = unique(.$chain)) %>% 
  ggplot() +
  aes(x = as.factor(generation), y = mean, group = chain, alpha = highlight == chain) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.1) +
  geom_line() + # geom_point() +
  facet_wrap(~highlight) +
  scale_alpha_discrete(range = c(.2, 1)) +
  theme(legend.position="none") +
  labs(x = "Generation", y = "Identifiability")
ident2

#--------------------------------------------------------
#' # Identifiability (models)
#--------------------------------------------------------

#' Does identifiability increase?

# full
mm.ident <- lmer(ident ~ generation + (1|chain), d, REML = F); summary(mm.ident)
# null
mm.ident.null <- lmer(ident ~ 1 + (1|chain), d, REML = F)
# comparison
anova(mm.ident.null, mm.ident)

#' Only with human data?

# full
mm.ident.h <- lmer(ident ~ generation + (1|chain), dh, REML = F); summary(mm.ident.h)
# null
mm.ident.null.h <- lmer(ident ~ 1 + (1|chain), dh, REML = F)
# comparison
anova(mm.ident.null.h, mm.ident.h)

#--------------------------------------------------------
#' # interesting ngrams (1)
#--------------------------------------------------------

#' frequency of ngrams per subject,
#' n of sequences containing an ngram

# all possible ngrams of size 2,3,4.
# 4 numbers representing the 4 syllables (ban, bi, ta, tin)
keys <- as.character(1:4)
# all possible combinations, including the boundary symbol (.)
pgrams <- combis(2, 4, c(keys, "."))

# combine all possible ngrams with all subjects
pgrams.subjs <- crossing(d$subjID, pgrams$ngram) %>% 
  rename(subjID = `d$subjID`, ngram = `pgrams$ngram`)

# combine all possible ngrams with all sequences.
# filter only sequences which contain the ngram in question.
# count the number of sequences per subject containing the given ngram.
# mark all missing ngrams with a frequency of 0
pgram.f <- d %>%
  mutate(seqOUT.boundaries = paste0(".", seqOUT, ".")) %>%
  select(subjID, seqOUT.boundaries) %>% 
  left_join(pgrams.subjs) %>% 
  filter(str_detect(seqOUT.boundaries, fixed(ngram))) %>% 
  count(subjID, ngram) %>% 
  right_join(pgrams.subjs) %>% 
  left_join(pgrams) %>% 
  mutate(n = ifelse(is.na(n), as.integer(0), n))

#--------------------------------------------------------
#' # interesting ngrams (2)
#--------------------------------------------------------

#' freq of ngrams in a corpus of 1,000,000 random sequences
#' generated using the same mehtod as for initial generations

# function generates shuffled sequences given one base sequence of symbs
stimuli <- function(symbs, len, trialn, seed=4) {
  # basic sequence to shuffle
  pool <- rep(symbs, len)
  # set seed for reproductibility
  set.seed(seed)
  # create shuffled sequences
  stim.seq <- character(trialn)
  for (i in seq(trialn)) {
    stim.seq[i] <- paste(sample(pool, length(pool), replace = F), collapse = "")
  }
  return(stim.seq)
}

#' running the shuffling function for 1M sequences takes ~18.643 seconds.
#' an alternative creating a tibble, then using purrr::map takes longer (~60s)
stims <- stimuli(keys, 3, 1000000)

#' add boundary symbols
stims <- paste(".", stims, ".", sep = "")

# saveRDS(stims, "hch2_baseline.rds")
# stims <- readRDS("hch2_baseline.rds")

#' placeholder for ngram frequencies in 1M-baseline
init.like <- pgrams %>%
  mutate(freq = NA_integer_)
#' run loop to find n of sequences in baseline containing a given ngram
#' takes ~90s
for (i in seq(nrow(pgrams))) {
  init.like$freq[i] <- length(grep(init.like$ngram[i], stims, fixed = T))
}

base.prob <- init.like %>% 
  mutate(prob = freq / length(stims))

# saveRDS(base.prob, "hch2_baseprob.rds")
# base.prob <- readRDS("hch2_baseprob.rds")

#--------------------------------------------------------
#' # interesting ngrams (3)
#--------------------------------------------------------

#' frequency ratio:
#' subjects' ngram frequency / relative to / 1M-baseline frequency

# base for all chains
pgram.base2 <- base.prob %>% 
  mutate(n = freq) %>% 
  select(gramsize, ngram, n)

# smoothing parameter when using the 10000000 baseline
minbase2 <- pgram.base2 %>% filter(n != 0) %>% select(n) %>% min() / 1000000
minbase2 <- minbase2*30

# ratio of ngram freq for each subject with respect to baseline
pgram.rat <- pgram.f %>% 
  left_join(pgram.base2, by = c("gramsize", "ngram")) %>%   # merge with simulated baseline
  mutate(n.subj = n.x) %>% 
  rename(f.subj = n.x, f.base2 = n.y) %>%
  filter(f.base2 != 0 | f.subj != 0) %>%                    # only keep ngrams present in subjects or simulated baseline
  left_join(pgram.base1, by = c("gramsize", "ngram", "chain")) %>%  # merge with 'real' baseline, = generation 0
  rename(f.base1.smooth = n) %>%
  mutate(f.subj.smooth = f.subj) %>% 
  mutate_each(funs((. + 1) / (max(.) + 1)), c(f.subj.smooth, f.base1.smooth)) %>%    # frequency with Laplace smoothing
  mutate_each(funs((. + minbase2) / (max(.) + minbase2)), c(f.base2, f.subj)) %>%    # frequency without smoothing
  mutate(logratio1 = log2(f.subj.smooth / f.base1.smooth),
         logratio2 = log2(f.subj / f.base2))

#--------------------------------------------------------
#' roadmap
#--------------------------------------------------------

# *** continue below
#' track:
#' gen.rat.lm >
#' pall.rat >
#' pgram.rat + pfeat.rat
#' pfeat.rat > pfeat.subj > pfeat.freq > p[onset|nucl|coda].subj > p[] + base2
#' pgram.rat > pgram.f > pgram.base2 > base.prob

#--------------------------------------------------------
#' over/under-represented ngrams
#--------------------------------------------------------

# place holder
# pgram.freq <- tbl_df(data.frame(subjID = character(0),
#                                 ngram = character(0),
#                                 n = integer(0), stringsAsFactors = F))
# system.time(
# for (i in pgrams$ngram) {
#   pgram.freq <- d %>%
#     group_by(subjID) %>%
#     filter(str_detect(outseq.b, fixed(i))) %>%
#     count(subjID) %>%
#     mutate(ngram = i) %>%
#     bind_rows(pgram.freq)
# }
# )
# saveRDS(pgram.freq, "../datuak/hch11/hch11_pgram.freq.rds")

# ngram frequency for each subject
# pgram.f <- tbl_df(pgrams[rep(seq(nrow(pgrams)), length(unique(d$subjID))),]) %>%
#   mutate(subjID = rep(unique(d$subjID), each = nrow(pgrams))) %>%
#   left_join(pgram.freq, by = c("ngram", "subjID")) %>%
#   mutate(n = replace(n, is.na(n), 0)) %>%
#   separate(col = subjID, into = c("chain", "generation"), "\\.", remove = F) %>%
#   mutate(generation = as.integer(generation))

# saveRDS(pgram.f, "../datuak/hch11/hch11_pgram.f.rds")
pgram.f <- readRDS("../datuak/hch11/hch11_pgram.f.rds")

# create 2 different frequency baselines to use as comparison:
# 1: the actual freq in generation 0 of each chain,
# 2: the freq in a corpus of 1,000,000 random sequences generated by the same mehtod

# base for all chains
pgram.base2 <- base.prob %>% 
  mutate(n = freq) %>% 
  select(gramsize, ngram, n)

# chain-specific base
pgram.base1 <- pgram.f %>% 
  filter(generation == 0) %>% 
  select(gramsize, ngram, n, chain)

# smoothing parameter when using the 10000000 baseline
minbase2 <- pgram.base2 %>% filter(n != 0) %>% select(n) %>% min() / 1000000
minbase2 <- minbase2*30

# View(plotgrams(as.character(ch), "^\\.", cutpoint = ctp, baseline = baseline)$data)

# ratio of ngram freq for each subject with respect to baseline
pgram.rat <- pgram.f %>% 
  left_join(pgram.base2, by = c("gramsize", "ngram")) %>%   # merge with simulated baseline
  mutate(n.subj = n.x) %>% 
  rename(f.subj = n.x, f.base2 = n.y) %>%
  filter(f.base2 != 0 | f.subj != 0) %>%                    # only keep ngrams present in subjects or simulated baseline
  left_join(pgram.base1, by = c("gramsize", "ngram", "chain")) %>%  # merge with 'real' baseline, = generation 0
  rename(f.base1.smooth = n) %>%
  mutate(f.subj.smooth = f.subj) %>% 
  mutate_each(funs((. + 1) / (max(.) + 1)), c(f.subj.smooth, f.base1.smooth)) %>%    # frequency with Laplace smoothing
  mutate_each(funs((. + minbase2) / (max(.) + minbase2)), c(f.base2, f.subj)) %>%    # frequency without smoothing
  mutate(logratio1 = log2(f.subj.smooth / f.base1.smooth),
         logratio2 = log2(f.subj / f.base2))

dfile <- paste("../datuak/hch11/hch11_", phase, "_pgram.rat.rds", sep = ""); dfile
saveRDS(pgram.rat, dfile)
