##################################################
# Tiny Overlap Finder
# based on the Minimal Generalisation Learner
# except I can't replicate the MGL output
# v0.1
##################################################

# -- fun -- #

# Hungarian orthography: replace characters in digraphs with their IPA equivalents or vice versa
transcribeIPA = function(string, direction){
  if (direction == 'single'){
    stringr::str_replace_all(string, c(
      'ccs' = 'cscs', 'ssz' = 'szsz', 'zzs' = 'zszs', 'tty' = 'tyty', 'ggy' = 'gygy', 'nny' = 'nyny', 'lly' = 'jj', 'cs' = 'č', 'sz' = 'ß', 'zs' = 'ž', 'ty' = 'ṯ', 'gy' = 'ḏ', 'ny' = 'ṉ', 'ly' = 'j', 's' = 'š', 'ß' = 's'))
  } else if (direction == 'double'){
    stringr::str_replace_all(string, c('s' = 'ß', 'š' = 's', 'ṉ' = 'ny', 'ḏ' = 'gy', 'ṯ' = 'ty', 'ž' = 'zs', 'ß' = 'sz', 'č' = 'cs'))
  }
}

# turn training set into mgl-compatible input
# legacy!
pairsIntoT = function(dat){
  dat |> 
  dplyr::mutate(
    # we need a binary category. if the log odds of the pair is above the mean, it's 'high', otherwise 'low'
    category = ifelse(log_odds >= mean(log_odds), 'high', 'low'),
    # transcribe input
    input = transcribeIPA(base, 'single'),
    # transcribe var1 if high and var2 if low
    output = dplyr::case_when(
      category == 'high' ~ transcribeIPA(form_1, 'single'),
      category == 'low' ~ transcribeIPA(form_2, 'single')
    ),
    # add TOF boundary markers
    output = glue::glue('#{output}#'),
    input = glue::glue('#{input}#')
  ) |> 
  dplyr::select(input,output)
}

# take two strings, find overlap from left, return overlap
# example: "abba" and "abbc" returns "abb"
# something and nothing returns nothing
str_overlap_left = function(str1, str2) {
  if (nchar(str1) == 0){
    return('')
  } else if (nchar(str2) == 0){
    return('')
  }
  else if (str1 == str2){
    return(str1)
  } else {
  # turn strings into vectors
  v1 = str_split_1(str1, '')
  v2 = str_split_1(str2, '')
  # take max length of shorter strings (overlap won't be more than that)
  max_length = min(length(v1), length(v2))
  # truncate strings
  v1 = v1[1:max_length]
  v2 = v2[1:max_length]
  # loop through strings, if they match, add to overlap, if not, break
  overlap = as.list(NULL)
  for (i in 1:max_length){
    if (v1[i] == v2[i]){
      overlap[[i]] = v1[i]
    } else {
      break
    }
  }
  # turn overlap back into string
  overlap = stringr::str_c(overlap, collapse = '')
  return(overlap)
  }
}

# take two strings, find overlap from right, return overlap
# example: "abba" and "cbba" returns "bba"
# something and nothing returns nothing
str_overlap_right = function(a, b){
  # reverse strings
  a2 = stringi::stri_reverse(a)
  b2 = stringi::stri_reverse(b)
  # find left overlap
  rev_right = str_overlap_left(a2, b2)
  # reverse results ;)
  stringi::stri_reverse(rev_right)
}

# delete part b of string a
str_delete = function(a, b){
  if (nchar(b) == 0){
    return(a)
  } else {
    stringr::str_replace(a, b, '') 
  }
}

# take training input, find rules: A -> B / C _ D
# you need ABCD across input-output pairs
formatTraining = function(dat){
  dat |>
  dplyr::rowwise() |>
  dplyr::mutate(
    # set up TOF boundary markers
    input = glue::glue('#{input}#'),
    output = glue::glue('#{output}#'),
    # find left overlap (C)
    c = str_overlap_left(input, output),
    # find what remains after C
    rem1 = str_delete(input, c),
    rem2 = str_delete(output, c),
    # find right overlap (D)
    d = str_overlap_right(rem1, rem2),
    # find what remains from input (A)
    a = str_delete(rem1, d),
    # and output (B)
    b = str_delete(rem2, d),
  ) |> 
  dplyr::select(-rem1,-rem2)
}

# cross rules on output of formatTraining
# now you have rules like
# x -> y / abc _ d
# and
# x -> y / fgc _ d
# but you want 
# x -> y / c _ d
# so you cross EVERY Input-Output pair with every other I-O pair and find overlaps in A,B,C,D

findLargerRules = function(t){
  # def left-hand side of cross
  r1 = t |> 
    dplyr::select(a,b,c,d) |> 
    dplyr::rename_all(~ glue::glue('{.}1'))
  
  # def right-hand side of cross
  r2 = t |> 
    dplyr::select(a,b,c,d) |> 
    dplyr::rename_all(~ glue::glue('{.}2'))
  
  # cross r with itself, mark rule pairs that do the same thing but have different contexts
  # where abc_d and fgc_d are different
  all_rules = tidyr::crossing(r1,r2) |> 
    dplyr::filter(a1 == a2 & b1 == b2) |> 
    dplyr::filter(!(c1 == c2 & d1 == d2)) |> 
    # string operations can do vectors, we don't want them to do that, so we do this rowwise
    dplyr::rowwise() |> 
    dplyr::mutate(
      # find the tiny overlap in the two C-s
      C = str_overlap_right(c1, c2),
      # and the two D-s
      D = str_overlap_left(d1, d2)
    ) |>
    dplyr::rename(
      # the rule's A and B remain the same
      A = a1, 
      B = b1
    ) |> 
    # take ABCD, drop duplicates
    dplyr::distinct(A,B,C,D) |> 
    dplyr::mutate(
      # set up rule
      rule = glue::glue('{A} → {B} / {C}__{D}')
    )
  
  return(all_rules)
}

# take output of formatTraining and output of findLargerRules and make rules and words set
# that is, take the now more general rules and join them back with the input-output pairs of words
getRulesAndWords = function(t,all_rules){
  
  # we cross all training I-O pairs with all rules
  rules_and_words = tidyr::crossing(t,all_rules) |> 
    dplyr::mutate(
      # we keep I-O pair if A and B match (so the rule applies to the pair)
      # and the rule's context overlaps with the input
      # so like the input is 'abc_d' and the rule is 'c_d' or 'bc_d' but not 'cba_d'
      # if CAD in the input is matches rule's context, the word is in the rule's scope
      scope = a == A & stringr::str_detect(c, glue::glue('{C}$')) & stringr::str_detect(d, glue::glue('^{D}')),
      # if the output matches the rule's output, the rule is a hit
      hit = scope & b == B
    ) |> 
    dplyr::filter(scope)
}

# take rules and words set and dplyr::count matches and scope for rules, then calc confidence and so on
# so now we know, for each rule, the n hits and n scope. we use these to calc rule metrics
getRuleStats = function(rules_and_words,alpha_upper,alpha_lower){
  
  # for each rule, scope is hits and misses, reliability is hits / scope
  all_rules_stats = rules_and_words |> 
    dplyr::count(rule, A, B, C, D, hit) |> 
    tidyr::pivot_wider(names_from = hit, values_from = n, values_fill = 0) |>
    dplyr::rename(hits = `TRUE`, misses = `FALSE`) |> 
    dplyr::mutate(
      scope = hits + misses,
      reliability = hits / scope
    )
  
  # now we go through each rule and
  all_rules_stats = all_rules_stats |>
    dplyr::rowwise() |>
    dplyr::mutate(
      # adjust reliability following Albright and Hayes and Mikheev
      adjusted_reliability = (hits + .5) / (scope + 1),
      # calculate the lower confidence value using the alpha lower and the t distribution.
      # Albright and Hayes are using the normal distribution. Mikheev uses t.
      t_lower = qt(p = (1 - alpha_lower)/2, df = scope - 1, lower.tail = F),
      # calculate the upper confidence value using the alpha upper and the t distribution
      t_upper = qt(p = (1 - alpha_upper)/2, df = scope - 1, lower.tail = F),
      # calculate the sample variance for the rule
      sample_variance = sqrt((adjusted_reliability * (1 - adjusted_reliability)) / scope),
      # get lower and upper conf limit using these
      lower_confidence_limit = adjusted_reliability - t_lower * sample_variance,
      upper_confidence_limit = adjusted_reliability + t_upper * sample_variance
    )
  
  return(all_rules_stats)
}

# find overlap between two contexts, c1 and c2
# the most general rule applies to ALL contexts such that the description is a -> b / _
# this means that "nothing" in fact is the superset of every other context!
compareContextsC = function(c1,c2){
  if(nchar(c1) == 0){
    return(T)
  } else {
    return(stringr::str_detect(c1, glue::glue('{c2}$')))
  }
}

# take rule1, rule2, and the original rules and words table
# for superset rule1 and its subset rule2 find residue reliability of rule1 for positions where rule2 doesn't apply
# basically we go back to the rules and words table and then use it to find the setdiff of rule1 and rule2 and count rows
getResidueReliability = function(rule1,rule2,rules_and_words){
  
  # this is the larger rule's scope
  r1_scope = rules_and_words |> 
    dplyr::filter(rule == rule1) |> 
    dplyr::select(input,output)
  
  # this is the larger rule's hits
  r1_hits = rules_and_words |> 
    dplyr::filter(
      rule == rule1, 
      hit
    ) |> 
    dplyr::select(input,output)
  
  # this is the smaller rule's scope
  r2_scope = rules_and_words |> 
    dplyr::filter(rule == rule2) |> 
    dplyr::select(input,output)
  
  # this is the smaller rule's hits
  r2_hits = rules_and_words |>
    dplyr::filter(
      rule == rule2, 
      hit
    ) |> 
    dplyr::select(input,output)
  
  # calc scope and hits for residue
  scope_residue = dplyr::anti_join(r1_scope,r2_scope, by = join_by(input, output))
  hits_residue = dplyr::anti_join(r1_hits,r2_hits, by = join_by(input, output))
  reliability_residue = (nrow(hits_residue) + .5) / (nrow(scope_residue) + 1)
  return(reliability_residue)
}

# impugn rules. 
# if rule1 applies to a superset context of rule2, maybe rule2 does all the work!
# if that's true, we want to penalise rule1.
# we figure that out by checking the reliability of the "residue rule" that is rule1 applying to contexts that rule2 doesn't apply to. read this one more time.
# if this reliability is pathetic then rule2 is doing all the work.
# this is a form of residualisation
impugnRules = function(all_rules_stats,rules_and_words){
  
  # build context
  all_rules_stats = all_rules_stats |>
    dplyr::mutate(
      context = glue::glue('{C}__{D}')
    )
  
  # we compare every rule to every rule, again!
  all_rules_stats_2 = all_rules_stats |> 
    dplyr::rename_all(~ glue::glue('{.}_2'))
  
  # we cross rules with rules
  rule_cross = tidyr::crossing(all_rules_stats,all_rules_stats_2) |>
    dplyr::filter(
      # do the two rules do the same thing?
      A == A_2, # same thing turns into...
      B == B_2, # same other thing,
      rule != rule_2, # but it's not the same rule!
      nchar(context) == min(nchar(context)) | nchar(context) < nchar(context_2), # we want the broader rules on the left. so we keep smaller string contexts (smaller context = applies to more forms) or the "everything context"*
      nchar(context_2) != min(nchar(context_2))
    ) |> 
    # we now have the rules in superset-subset relationships lined up next to each other
    # we go through the rule pairs and compare contexts
    dplyr::rowwise() |> 
    dplyr::mutate(
      keep = compareContextsC(C_2, C)
    ) |> 
    dplyr::filter(keep)
  
  # *if like a rule 1 is a -> b / cd_e and rule 2 is a -> b / d_e then rule 2 will apply everywhere where rule 1 applies, since words ending in 'cd' is a subset of words ending in 'd'
  
  # rule 1 is always the bigger rule, since context 'a$' will always superset context 'ba$'
  # now I need to go back to rb and calc the rel(C1-C2) = (hits(C1) - hits(C2))/ (scope(C1) - scope(C2))
  # for that I need to grab hits and scope words for C1 and C2 from all rules w/ words
  
  # now we calc residue reliability
  rule_cross = rule_cross |>
    dplyr::mutate(
      adjusted_residue_reliability = purrr::map2_dbl(rule, rule_2, ~ getResidueReliability(.x,.y,rules_and_words))
    )
  
  # the residue rule also has conf limits (where we imagine that it applies to infinite n forms)
  big_rules = rule_cross |> 
    dplyr::mutate(
      # we calc the UPPER confidence limit for the residue rule
      upper_residue_confidence_limit = adjusted_residue_reliability + t_upper * sample_variance,
      # if the LOWER confidence limit for the superset rule is LOWER than the UPPER confidence limit for the residue rule, then it's bad and we keep the residue rule, since that's all there is to the superset rule, really.
      # if the LOWER confidence limit for the superset rule is HIGHER than the UPPER confidence limit for the residue rule, then we keep the superset rule, since it's doing work.
      # we compare UPPER with LOWER to get a conservative estimate
      # this is all Albright and Hayes of course
      impugned_lower_confidence_limit = ifelse(
        upper_residue_confidence_limit < lower_confidence_limit, upper_residue_confidence_limit, upper_confidence_limit)
    ) |> 
    # all we get from this is the impugned_Lower_conf_limit, which is the final confidence of the rule
    dplyr::distinct(A,B,C,D,rule,scope,hits,reliability,adjusted_reliability,lower_confidence_limit,impugned_lower_confidence_limit)
  
  # we now have the big rules. we take the original rules, drop the big rules, add an upper_res_conf_lim (which, for these, will be the same as the lower_conf) and then bind rows
  # we keep both the bigger rule and the smaller rule, we don't prune rules
  new_rules = all_rules_stats |>
  dplyr::filter(!rule %in% big_rules$rule) |>
  dplyr::mutate(impugned_lower_confidence_limit = lower_confidence_limit) |>
  dplyr::select(A,B,C,D,rule,scope,hits,reliability,adjusted_reliability,lower_confidence_limit,impugned_lower_confidence_limit) |>
  dplyr::bind_rows(big_rules)

  return(new_rules)
}

# add examples and exceptions
# we now have rules and rule metrics. we want to know which words are in the rule's scope and which words the rule applies to.
# we go back to the words and rules table and fish out the words for each rule
addWordsToRules = function(rules,rules_and_words){
  
  # take the orthographic forms and the hit col for each rule (so we know if the rule applies or not)
  rules_w_e = rules_and_words |>
    dplyr::select(orth,hit,rule) |>
    dplyr::group_by(rule,hit) |>
    # we glue together the hits and the exceptions into two strings
    tidyr::nest() |>
    dplyr::mutate(
        cat = dplyr::case_when(
          hit ~ 'examples',
          !hit ~ 'exceptions'
          ),
        words = purrr::map_chr(data, ~ .x |> dplyr::pull(orth) |> paste(collapse = ', '))
    ) |>
    dplyr::ungroup() |>
    # we get the rule and the words
    dplyr::select(rule,cat,words) |>
    # we put the hits and exceptions next to each other
    tidyr::pivot_wider(names_from = cat, values_from = words, values_fill = '') |>
    dplyr::select(rule,examples,exceptions)
  
  # we merge the rules which have all the metrics and the rules with the hits and exception lists
  left_join(rules,rules_w_e, by = join_by(rule))
  # tadaah
}

# -- run the whole thing -- #

# the ahem tiny overlap finder needs a training set and two alpha values
# I interpersed the code with quotes from the original papers.
# it's all Albright and Hayes except the bit where they refer the reader to Mikheev
tof = function(training,alpha_lower,alpha_upper){
  
  # check for input column in training
  if(!'input' %in% colnames(training)){
    stop('training data must have a column named "input"')
  }
  # check for output column in training
  if(!'output' %in% colnames(training)){
    stop('training data must have a column named "output"')
  }
  # check for alpha_lower
  if(!is.numeric(alpha_lower)){
    stop('you have to define alpha_lower as a number between 0-1')
  }
  # check for alpha_upper
  if(!is.numeric(alpha_upper)){
    stop('you have to define alpha_upper as a number between 0-1')
  }
  
  # Modeling English Past Tense Intuitions with Minimal Generalization Albright and Hayes 2002
  # also: Automatic Rule Induction for Unknown-Word Guessing Mikheev 1997
  
  # "Formally, the structural change can be represented in the format A → B, and the context in the format / C__D, to yield word-specific rules like those in (2). (The symbol ‘#’ stands for a word boundary.)"
  
  # "The exact procedure for finding a word-specific rule is as follows: given an input pair (X, Y), the model first finds the maximal left-side substring shared by the two forms (e.g., #ms), to create the C term (left side context). The model then exam- ines the remaining material and finds the maximal substring shared on the right side, to create the D term (right side context). The remaining material is the change; the non-shared string from the first form is the A term, and from the second form is the B term."
  
  d = formatTraining(training)
  
  # [Now you need to find overlaps between structural descriptions of rules by comparing every rule to every other rule]
  
  d2 = findLargerRules(d)
  
  # combine rules with words
  
  d3 = getRulesAndWords(d,d2)
  
  # [Now you find scope, and hits for each rule and then calculate lower conf following Mikhailev 1997]
  
  # First, we determine the number of forms in the training data that meet the structural description of the rule (for A → B / C__D, these are the forms that contain CAD). This number is the scope of the rule. The hits of the rule is the number of forms that it actually derives correctly. The reliability of a rule is simply the ratio of its hits to its scope.
  
  # Mikhaleev:
  # The estimate is a good indicator of the rule accuracy but it frequently suffers from large estimation error due to insufficient training data. For example, if a rule was found to apply just once and the total number of observations was also one, its estimate p has the maximal value (1) but clearly this is not a very reliable estimate. We tackle this problem by calculating the lower confidence limit 71"L for the rule estimate, which can be seen as the minimal expected value of/~ for the rule if we were to draw a large number of samples. 
  
  d4 = getRuleStats(d3, alpha_upper, alpha_lower)
  
  # [ we skip phonological rules. for this particular set of uh datasets I've done morpho-phonology beforehand (like, I've removed VH from -ik, split others into suffix classes)]
  
  # formalizing this intuition, we propose a re- finement of the way that confidence is calculated, in order to diagnose when a subpart of a generali- zation is doing most of the work of the larger gen- eralization. When we consider the confidence of a context C associated with a change A → B, we must consider every other context C′ associated with A → B, checking to see whether C′ covers a subset of the cases that C covers. In the present case, when we assess the confidence of adding [-t] after any consonant, we would check all of the other rules adding [-t], including the one that adds [-t] after obstruents. For each C′ that covers a sub- set of C,we must ask whether the rule A→B/C′ is actually “doing most of the work” of the larger rule A→B/C.
  
  # To find out if the smaller rule is doing most of the work, we calculate how well the larger rule (C) performs outside the area covered by the smaller rule (C′). The reliability of the residue area (C -C') is calculated as follows: rel(C-C') = (hits(C) - hits(C'))/ (scope(C) - scope(C')).
  
  # From the reliability of this residue area (C – C′), we can then calculate its confidence, using confidence limit statistics in a way similar to that described above in section 3.2. However, there is a crucial difference: when we are assessing whether a rule explains enough cases to be trustable, we are interested in the denseness of cases within the gen- eralization. But when we are assessing whether a rule offers an improvement over a subpart, we are interested in the sparseness of cases in the residue outside of the subpart. Therefore, when calculating the confidence of the residue, we must use the up- per confidence limit rather than the lower confi-dence limit.
  
  # If the upper confidence limit of the reliability of the residue (C – C′) is lower than the lower confi- dence limit of the reliability of the larger context (C), then we can infer that the smaller rule (A→B/C′) is doing most of the work of the larger rule (A → B / C). Therefore, we penalize the larger rule by replacing its confidence value (Lower confidence(C)) with the confidence value of the residue (Upper confidence(C – C′)). We call this penalty impugnment, because the validity of the larger rule is being called into question by the smaller rule. Impugnment is carried out for all contexts of all rules.
  
  d5 = impugnRules(d4, d3)
  
  # add examples and exceptions to each rule
  
  d6 = addWordsToRules(d5, d3)
  
  return(d6)
}
