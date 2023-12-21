# preprocessing for the three variables

setwd('~/Github/packages/worddist')

library(tidyverse)
library(lme4)
library(hunspell)

source('code/helper.R')

# -- read -- #

training_ik = read_tsv('~/Github/Racz2024/resource/real_words/ik_verbs/ikes_pairs_webcorpus2.tsv')
training_ep = read_tsv('~/Github/Racz2024/resource/real_words/epenthetic_stems/epenthesis_pairs_webcorpus2.tsv')
training_vh = read_tsv('~/Github/Racz2024/resource/real_words/front_harmony/fh_pairs_webcorpus2.tsv')
test = read_tsv('~/Github/Racz2024/exp_data/baseline/baseline_tidy_proc.tsv')

# -- test -- #

test2 = test %>% 
  mutate(
    string = case_when(
      variation == 'lakok/lakom' ~ base %>% 
        str_replace('ik$', '') %>% 
        transcribeIPA('single'),
      variation == 'cselekszenek/cselekednek' ~ base %>% 
        transcribeIPA('single') %>%
        str_extract('.*(?= \\/)') %>% 
        str_replace('ik', ''),
      variation == 'hotelban/hotelben' ~ base %>% 
        transcribeIPA('single')
    ),
    nchar = nchar(string)
  )

# -- training -- #

# - ik - #

training_ik2 = training_ik  %>% 
  mutate(
    variation = 'lakok/lakom', # this is arbitrarily wrong, the order of k and m forms (k is form1 m is form2) is the same as in test
    string = base %>% 
      str_replace('ik$', '') %>% 
      transcribeIPA('single'),
    nchar = nchar(string)
  ) %>% 
  filter(
    log(lemma_freq_corrected) > quantile(log(lemma_freq_corrected), .25),
    freq_1 > 1,
    freq_2 > 1,
    nchar < quantile(nchar, .9)
  ) %>% 
  select(base,variation,lemma_freq_corrected,string,log_odds)

# - ep - #

fit1 = glmer(cbind(freq_1,freq_2) ~ 1 + (1|base) + (1|xpostag), family = binomial, data = training_ep)

training_ep2 = training_ep %>% 
  mutate(
    string = base %>% 
      str_replace('ik$', '') %>% 
      transcribeIPA('single'),
    nchar = nchar(string)
  ) %>% 
  distinct(base,variation,lemma_freq_corrected,string)

training_ep3 = as.data.frame(ranef(fit1)$base) %>% 
  rownames_to_column() %>% 
  rename(base = rowname, ranef = `(Intercept)`) %>% 
  inner_join(training_ep2)

# -- vh -- #

# adapted from RaczRebrus2024
training_vh2 = training_vh %>% 
  filter(
    str_detect(base, '[xw]', negate = TRUE),
    hunspell_check(base, dict = dictionary('hu_HU')),
    str_detect(base, '(átles|enyv$|amely|hárem|moment|perc$|szeg$|nem$|blues$|assembly$|jegy$|szer$|cosec$|tett$|keksz$|dressz$|mez$|szenny$|szex$|szerv$|terv$|elv$|est$|ent$|elt$|sejt$|les$|csekk$|szleng$|kedd$|szerb$|szesz$|hely$|fej$|jel$|rend$|szocdem$|nyelv$|kert$|test$|teszt$|szett$|csepp$|perm$|mell$|csel$|kel$|szem$|kedv$|hegy$|szent$|meccs$|vers$|meggy$|borzeb|sosem|nedv$|necc$|tej$|segg$|csend$|seb$|kommersz|kokett|móres|goeb|gazeb)', negate = T)
  ) %>% 
  mutate(
    base = str_replace(base, 'boyler', 'bojler')
  )

fit2 = glmer(cbind(freq_1,freq_2) ~ 1 + (1|base) + (1|xpostag), family = binomial, data = training_vh2)

training_vh2 = training_vh2 %>% 
  mutate(
    string = base %>% 
      transcribeIPA('single'),
    nchar = nchar(string)
  ) %>% 
  distinct(base,variation,lemma_freq_corrected,string)

training_vh3 = as.data.frame(ranef(fit2)$base) %>% 
  rownames_to_column() %>% 
  rename(base = rowname, ranef = `(Intercept)`) %>% 
  inner_join(training_vh2)

# -- combine -- #

training2 = bind_rows(
  training_ik2,
  training_ep3,
  training_vh3
)

# -- write -- #

write_tsv(test2, 'dat/training_sets/test_set.tsv')
write_tsv(training2, 'dat/training_sets/training_set.tsv')
