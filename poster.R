setwd('~/Github/Racz2024b')

set.seed(1338)
options(knitr.kable.NA = "")

library(tidyverse)
library(glue)
library(magrittr)
library(knitr)
library(ggthemes)
library(lme4)
library(patchwork)
library(ggrepel)
library(stringdist)

source('code/helper.R')

tests = read_tsv('dat/tests.tz')
# word_distance = read_tsv('dat/alignments/word_distances.tsv')

# phonology/morphology

p1 = tests |>  
  filter(str_detect(variation, '(lakok|cselek)')) |> 
  HRV() |> 
  mutate(
    derivational = str_replace(derivational,'sz','s') |> 
      fct_relevel('-sik','-zik'),
    `log odds` = log((resp1+1)/(resp2+1))
  ) |> 
  ggplot(aes(derivational, `log odds`, colour = as.factor(nsyl))) +
  gghalves::geom_half_violin(side = 'r') +
  gghalves::geom_half_boxplot(side = 'r', width = .15, position = position_dodge(width = .75)) +
  gghalves::geom_half_point(side = 'l') +
  theme_bw() +
  scale_colour_grey() +
  xlab('derivational suffix') +
  coord_flip() +
  labs(colour = 'number of\nstem\nsyllables') +
  facet_wrap( ~ variation) +
  theme(strip.background = element_blank())

p2 = tests |>  
  filter(str_detect(variation, 'hotel')) |> 
  HRV() |> 
  mutate(
    `log odds` = log((resp1+1)/(resp2+1))
  ) |> 
  ggplot(aes(vowel, `log odds`)) +
  gghalves::geom_half_violin(side = 'r') +
  gghalves::geom_half_boxplot(side = 'r', width = .15, position = position_dodge(width = .75)) +
  gghalves::geom_half_point(side = 'l') +
  theme_bw() +
  scale_colour_grey() +
  xlab('final stem vowel') +
  coord_flip() +
  labs(colour = 'number of\nstem\nsyllables') +
  facet_wrap( ~ variation) +
  theme(strip.background = element_blank())

p2 + p1 + plot_layout(widths = c(1,2))
ggsave('figures/poster1.pdf', width = 12, height = 3)

# distributions

tests |>
  HRV() |> 
  mutate(
    `log odds` = log((resp1+1)/(resp2+1)) |> 
      scales::rescale()
  ) |> 
  select(base,variation,`log odds`,mgl,gcm,knn) |> 
  pivot_longer(-c(base,variation)) |> 
  mutate(
    name = fct_relevel(name, 'log odds', 'knn')
  ) |> 
  ggplot(aes(value)) +
  geom_histogram(alpha = .75) +
  facet_wrap( ~ variation + name) +
  theme_bw() +
  xlab('values (rescaled)') +
  theme(strip.background = element_blank())

ggsave('figures/poster2.pdf', width = 9, height = 6)

# results: best model fits
# this is the worst possible way of doing this.

broom::tidy(glm(cbind(resp1,resp2) ~ 1 + knn, data = tests[tests$variation == 'lakok/lakom',], family = binomial), conf.int = T) |> filter(term == 'knn') |> select(conf.low,conf.high)
broom::tidy(glm(cbind(resp1,resp2) ~ 1 + knn, data = tests[tests$variation == 'cselekszenek/cselekednek',], family = binomial), conf.int = T) |> filter(term == 'knn') |> select(conf.low,conf.high)
broom::tidy(glm(cbind(resp1,resp2) ~ 1 + knn, data = tests[tests$variation == 'hotelban/hotelben',], family = binomial), conf.int = T) |> filter(term == 'knn') |> select(conf.low,conf.high)

broom::tidy(glm(cbind(resp1,resp2) ~ 1 + gcm, data = tests[tests$variation == 'lakok/lakom',], family = binomial), conf.int = T) |> filter(term == 'gcm') |> select(conf.low,conf.high)
broom::tidy(glm(cbind(resp1,resp2) ~ 1 + gcm, data = tests[tests$variation == 'cselekszenek/cselekednek',], family = binomial), conf.int = T) |> filter(term == 'gcm') |> select(conf.low,conf.high)
broom::tidy(glm(cbind(resp1,resp2) ~ 1 + gcm, data = tests[tests$variation == 'hotelban/hotelben',], family = binomial), conf.int = T) |> filter(term == 'gcm') |> select(conf.low,conf.high)

broom::tidy(glm(cbind(resp1,resp2) ~ 1 + mgl, data = tests[tests$variation == 'lakok/lakom',], family = binomial), conf.int = T) |> filter(term == 'mgl') |> select(conf.low,conf.high)
broom::tidy(glm(cbind(resp1,resp2) ~ 1 + mgl, data = tests[tests$variation == 'cselekszenek/cselekednek',], family = binomial), conf.int = T) |> filter(term == 'mgl') |> select(conf.low,conf.high)
broom::tidy(glm(cbind(resp1,resp2) ~ 1 + mgl, data = tests[tests$variation == 'hotelban/hotelben',], family = binomial), conf.int = T) |> filter(term == 'mgl') |> select(conf.low,conf.high)

with(summary(glm(cbind(resp1,resp2) ~ 1 + knn, data = tests[tests$variation == 'lakok/lakom',], family = binomial)), 1 - deviance/null.deviance)
with(summary(glm(cbind(resp1,resp2) ~ 1 + gcm, data = tests[tests$variation == 'lakok/lakom',], family = binomial)), 1 - deviance/null.deviance)
with(summary(glm(cbind(resp1,resp2) ~ 1 + mgl, data = tests[tests$variation == 'lakok/lakom',], family = binomial)), 1 - deviance/null.deviance)

with(summary(glm(cbind(resp1,resp2) ~ 1 + knn, data = tests[tests$variation == 'cselekszenek/cselekednek',], family = binomial)), 1 - deviance/null.deviance)
with(summary(glm(cbind(resp1,resp2) ~ 1 + gcm, data = tests[tests$variation == 'cselekszenek/cselekednek',], family = binomial)), 1 - deviance/null.deviance)
with(summary(glm(cbind(resp1,resp2) ~ 1 + mgl, data = tests[tests$variation == 'cselekszenek/cselekednek',], family = binomial)), 1 - deviance/null.deviance)

with(summary(glm(cbind(resp1,resp2) ~ 1 + knn, data = tests[tests$variation == 'hotelban/hotelben',], family = binomial)), 1 - deviance/null.deviance)
with(summary(glm(cbind(resp1,resp2) ~ 1 + gcm, data = tests[tests$variation == 'hotelban/hotelben',], family = binomial)), 1 - deviance/null.deviance)
with(summary(glm(cbind(resp1,resp2) ~ 1 + mgl, data = tests[tests$variation == 'hotelban/hotelben',], family = binomial)), 1 - deviance/null.deviance)
