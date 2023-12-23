# preprocessing for the three variables

setwd('~/Github/Racz2024b')

library(tidyverse)
library(lme4)
library(hunspell)

source('code/helper.R')

# -- read -- #

training_ik = read_tsv('~/Github/Racz2024/resource/real_words/ik_verbs/ikes_pairs_webcorpus2.tsv')
training_ep = read_tsv('dat/training_sets/cselekszenek_alt_training_set.tsv')
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
  mutate(
    median_val = median(log_odds),
    category = case_when(
      log_odds >= median_val ~ 'high',
      log_odds < median_val ~ 'low'
    )
  ) %>% 
  select(base,variation,lemma_freq_corrected,string,category)
  
# - ep - #

igekotok = '^(be|ki|fel|föl|túl|le|meg|át|oda|vissza|el|külön|végig|hátra|hozzá|össze|rá|tovább|ide|oda|utána|közbe|szét)'

training_ep2 = training_ep %>% 
  filter(
    !lemma %in% c('kszik','latszik'),
    str_detect(lemma, '(ll|zz)ik$', negate = TRUE),
    str_detect(lemma, igekotok, negate = TRUE),
    lemma_freq > quantile(lemma_freq, .5),
    nchar(lemma) < 11 # so the truncated lemma is max as long as the longest string in test + 
         ) %>% 
  mutate(
  string = lemma %>% 
    str_replace('x', 'ksz') %>% 
    str_replace('w', 'v') %>% 
    transcribeIPA('single') %>% 
    str_replace('ik$', ''),
  category = case_when(
    stem_type == 'cc' ~ 'high',
    stem_type == 'cvc' ~ 'low'
  )
) %>% 
  mutate(
    variation = 'cselekszenek/cselekednek',
    lemma_freq_corrected = lemma_freq,
    base = lemma,
    category = case_when(
      str_detect(base, '[oóe][dz]ik$') ~ 'low',
      T ~ category # I have no clue why this ended up here.
    )
  ) %>% 
  select(base, variation, lemma_freq_corrected, string, category)

# double check
training_ep2 %>% 
  filter(category == 'high') %>% 
  pull(base) %>% 
  sort()

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
  inner_join(training_vh2) %>% 
  mutate(
    category = case_when(
      ranef > median(ranef) ~ 'high',
      ranef <= median(ranef) ~ 'low'
    )
  ) %>% 
  select(base,variation,lemma_freq_corrected,string,category)

# -- combine -- #

training2 = bind_rows(
  training_ik2,
  training_ep2,
  training_vh3
)

# -- write -- #

write_tsv(test2, 'dat/training_sets/test_set.tsv')
write_tsv(training2, 'dat/training_sets/training_set.tsv')

