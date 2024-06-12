
#' ---
#' title: "Calculate validation for 1500 circRNAs"
#' author: "Marieke Vromman"
#' output: 
#'    html_document:
#'       toc: TRUE
#'       toc_float: TRUE
#'       theme: paper
#'       df_print: paged
#'       highlight: tango
#' ---
#' 

#' # File set-up

#' ## Set working directory to current directory
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

#' ## Load standard libraries and resolves conflicts
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')


#' ## Read data
#' table with all measured Cq values
cq = read_tsv("../data/details/selected_circ_cq.txt")
cq

#' # Label each circ for three methods
#' ## qPCR validation

cq = cq %>%
  # qPCR_val is 'pass' for all circ that have at least one replicate > 10
  mutate(qPCR_val = ifelse((Cq_min_untreated < 10) & (Cq_max_untreated < 10), 'fail', 'pass'),
         qPCR_val = ifelse((is.na(Cq_min_untreated)) & (is.na(Cq_max_untreated)), 'fail', qPCR_val),
         # qPCR_val_detail describes the reason why a circ fails this validation step if it does
         qPCR_val_detail = ifelse((Cq_min_untreated < 10) & (Cq_max_untreated < 10), 'all_Cq_under_10', 'pass'),
         qPCR_val_detail = ifelse(is.na(Cq_min_untreated) & is.na(Cq_max_untreated), 'all_NA', qPCR_val_detail))

cq

#' statistics
cq %>% count(qPCR_val)
cq %>% count(qPCR_val_detail)


#' ## RR val

cq = cq %>%
  # calculate delta cq
  # first make temp 47 for NAs
  mutate(Cq_min_treated_tmp = Cq_min_treated,
         Cq_min_treated = ifelse(is.na(Cq_min_treated), 47, Cq_min_treated),
         Cq_diff = Cq_min_treated - Cq_max_untreated) %>% 
  mutate(# fail if diff > 3
         RR_val = ifelse(Cq_diff > 3, "fail", 'pass'),
         RR_val_detail = ifelse(RR_val == 'fail', "FP (Cq_diff > 3)", 'pass'),
         # when Cq untreated is high => put NA
         RR_val = ifelse(Cq_min_untreated > 32, NA, RR_val),
         RR_val_detail = ifelse(Cq_min_untreated > 32, 'out_of_range', RR_val_detail),
         # if qPCR already failed => RR also becomes 'fail'
         RR_val = ifelse(qPCR_val == "fail", 'fail', RR_val),
         RR_val_detail = ifelse(qPCR_val == "fail", "fail_qPCR", RR_val_detail)) %>%
  # clean up NAs for 47
  mutate(Cq_diff = ifelse(Cq_min_treated == 47, NA, Cq_diff),
         Cq_min_treated = Cq_min_treated_tmp) %>%
  select(-Cq_min_treated_tmp) %>%
  #' to clean up dataframe: also remove Cq_diff for 'all_Cq_under_10' group and out_of_range group
  mutate(Cq_diff = ifelse(qPCR_val_detail == 'all_Cq_under_10', NA, Cq_diff),
         Cq_diff = ifelse(RR_val_detail == 'out_of_range', NA, Cq_diff))

cq

#' statistics
cq %>% count(RR_val)
cq %>% count(RR_val_detail)


#' ## Amplicon sequencing
cq_amp = read_tsv('../data/details/amp_seq.txt')

#' statistics
cq_amp %>% count(amp_seq_val)
cq_amp %>% count(amp_seq_val_detail)


#' add to cq dataframe
cq = cq %>% 
  left_join(cq_amp %>% unique()) %>%
  mutate(amp_seq_val_detail = ifelse(is.na(amp_seq_val_detail), "not_included", amp_seq_val_detail))

#' statistics
cq %>% count(amp_seq_val)
cq %>% count(amp_seq_val_detail)


#' ## Add a compound validation 'decision' for each circ
cq = cq %>%
  mutate(RR_val_tmp = ifelse(RR_val_detail == 'out_of_range', 'fail', RR_val)) %>% # when circ cannot be measured for RR => count as fail
  mutate(compound_val = ifelse(qPCR_val == "pass" & RR_val_tmp == "pass" & amp_seq_val == "pass", "pass", "NA"),
         compound_val = ifelse(qPCR_val == 'fail' | RR_val_tmp == "fail" | amp_seq_val == 'fail', 'fail', compound_val),
         compound_val = ifelse(is.na(amp_seq_val), NA, compound_val)) %>% 
  select(-RR_val_tmp)

#' statistics
cq %>% count(compound_val)

#' ## Add circRNA annotation and save as Supplementary Table 3
all_circ = read_tsv("../data/Supplementary_Table_2_all_circRNAs.txt")

all_circ

cq = cq %>%
  left_join(all_circ) 

cq %>%
  write_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')

#' # Calculate val rates (precision) per method per tool
perc_val = cq %>% 
  select(tool, circ_id, cell_line, count_group, qPCR_val, RR_val, RR_val_detail, amp_seq_val, amp_seq_val_detail, compound_val) %>%
  group_by(tool, count_group) %>%
  summarise(nr_qPCR_total = n(),
            nr_qPCR_fail = sum(qPCR_val == 'fail'),
            nr_qPCR_val = sum(qPCR_val == 'pass'),
            nr_RR_total = n() - sum(is.na(RR_val)),  # here NA are the ones that have are 'out_of_range'
            nr_RR_fail = sum(RR_val == "fail", na.rm = T),
            nr_RR_val = sum(RR_val == "pass", na.rm = T),
            nr_amp_total = n() - sum(is.na(amp_seq_val)),  # here NA are the ones 'not_included'
            nr_amp_fail = sum(amp_seq_val == "fail", na.rm = T),
            nr_amp_val = sum(amp_seq_val == "pass", na.rm = T),
            nr_compound_total = n() - sum(is.na(amp_seq_val)), # here NA are the ones 'not_included'
            nr_compound_fail = sum(compound_val == "fail", na.rm = T),
            nr_compound_val = sum(compound_val == "pass", na.rm = T)) %>%
  mutate(perc_qPCR_val = nr_qPCR_val/nr_qPCR_total,
         perc_RR_val = nr_RR_val/nr_RR_total,
         perc_amp_val = nr_amp_val/nr_amp_total,
         perc_compound_val = nr_compound_val/nr_compound_total) %>%
  ungroup()

perc_val


#' ## Save dataframe (Sup Table 6A)

perc_val %>% write_tsv('../data/Supplementary_Table_6A_precision_values.txt')



#' # Calculate sensitivity
#' ## Based on validated set (n = 949)

#' first calculate the total nr of validated circ per count group

nr_val = cq %>% 
  # get set of uniquely validated circRNAs
  filter(compound_val == 'pass') %>%
  select(circ_id, cell_line, count_group_median) %>% unique() %>%
  count(count_group_median) %>%
  rename(nr_expected = n)

nr_val

#' then calculate sensitivity by dividing nr of circ found by total
sens = cq %>% 
  # get set of uniquely validated circRNAs
  filter(compound_val == 'pass') %>%
  select(circ_id, cell_line, count_group_median) %>% unique() %>%
  # check witch tools have detected these
  left_join(all_circ %>%
              select(tool, circ_id, cell_line) %>% unique()) %>%
  group_by(tool, count_group_median) %>% 
  summarise(nr_detected = n()) %>%
  left_join(nr_val) %>%
  mutate(sensitivity = nr_detected/nr_expected) %>%
  ungroup()

sens


#' ## Based on theoretical nr of TPs => extrapolated sensitivity
sens_2 = perc_val %>% 
  # use perc_compound_val
  select(tool, count_group, perc_compound_val) %>%
  # get the number of detected circRNAs for each cell line and tool, per count group
  left_join(all_circ %>%
              group_by(tool, count_group) %>%
              summarise(total_n = n())) %>%
  # calculate the theoretical nr of validated circ
  mutate(extrapolated_sensitivity = perc_compound_val * total_n) 

sens_2

#' take both types of sensitivity together and fix Sailfish-circ count group thing
sens = sens %>% 
  filter(!tool == 'Sailfish-cir') %>% 
  full_join(sens_2 %>% filter(!tool == 'Sailfish-cir'),
            by = c('count_group_median' = 'count_group', 'tool')) %>%
  bind_rows(sens %>% 
              filter(tool == 'Sailfish-cir') %>%
              full_join(sens_2 %>% filter(tool == 'Sailfish-cir') %>%
                        select(-count_group),
                        by = c('tool')))

sens

#' ## Save dataframe (Sup Table 6B)
sens %>% write_tsv('../data/Supplementary_Table_6B_sensitivity_values.txt')


#' # Make dataframe with all info per tool and ranking

val_sens_df = perc_val %>%
  select(tool, count_group, perc_qPCR_val, perc_RR_val, perc_amp_val, perc_compound_val) %>%
  full_join(sens %>% rename(count_group = count_group_median) %>%
              select(-nr_detected, -nr_expected, -total_n) %>%
              mutate(count_group = ifelse(tool == "Sailfish-cir", 'no_counts', count_group))) %>%
  full_join(all_circ %>%
              group_by(tool, count_group) %>% count()) %>%
  group_by(count_group) %>%
  mutate(qPCR_precision_rank = dense_rank(desc(perc_qPCR_val)),
         RR_precision_rank = dense_rank(desc(perc_RR_val)),
         amp_precision_rank = dense_rank(desc(perc_amp_val)),
         compound_precision_rank = dense_rank(desc(perc_compound_val)),
         sensitivity_rank = dense_rank(desc(sensitivity)),
         nr_circ_rank = dense_rank(desc(n)),
         extrapolated_sensitivity_rank = dense_rank(desc(extrapolated_sensitivity))) %>%
  ungroup() %>%
  filter(!is.na(n)) %>% 
  select(tool, count_group, n,	nr_circ_rank, perc_qPCR_val, qPCR_precision_rank,
         perc_RR_val, RR_precision_rank, perc_amp_val, amp_precision_rank,
         perc_compound_val, compound_precision_rank, sensitivity, sensitivity_rank,
         extrapolated_sensitivity, extrapolated_sensitivity_rank) %>%
  arrange(desc(count_group), tool) %>%
  rename(nr_circ_detected = n, qPCR_precision = perc_qPCR_val,
         RR_precision = perc_RR_val, amp_precision = perc_amp_val,
         compound_precision = perc_compound_val) %>% 
  mutate(across(c(qPCR_precision, RR_precision, amp_precision,
                  compound_precision, sensitivity), round, 4),
         extrapolated_sensitivity = round(extrapolated_sensitivity, 0))

val_sens_df

val_sens_df %>%
  write_tsv('../data/Supplementary_Table_6_tool_ranking.txt')


#' # RNase R enrichment calculations

#' calculate enrichment

all_circ_treated = read_tsv('../data/Supplementary_Table_4_all_circRNAs_treated.txt')

treatment = all_circ %>%
  rename(count_UT = BSJ_count) %>% 
  select(chr, start, end, strand, cell_line, tool, circ_id, circ_id_strand, count_UT, count_group) %>%
  full_join(all_circ_treated %>% 
              rename(count_T = BSJ_count) %>%
              select(chr, start, end, strand, cell_line, tool, circ_id, circ_id_strand, count_T))

treatment

treatment = treatment %>%
  left_join(read_tsv('../data/details/nr_reads.txt') %>% 
              select(-RNA_id) %>%
              pivot_wider(names_from = RNaseR, values_from = nr_reads) %>%
              rename(nr_reads_T = treated, nr_reads_UT = untreated)) %>%
  # calculate CPMs and enrichment factor (Sailfish-cir is already TPM)
  mutate(cpm_T = ifelse(tool == 'Sailfish-cir', count_T, 1000000 * count_T/nr_reads_T), 
         cpm_UT = ifelse(tool == 'Sailfish-cir', count_UT, 1000000 * count_UT/nr_reads_UT),
         enrichment = cpm_T/cpm_UT,
         enrichment_bin = ifelse(enrichment > 1, 'enriched', 'not enriched'),
         enrichment_bin = ifelse(is.na(count_T), 'count treated NA', enrichment_bin),
         enrichment_bin = ifelse(is.na(count_UT), 'count untreated NA', enrichment_bin))

treatment %>% write_tsv('../data/Supplementary_Table_5_RNase_R_enrichment_seq.txt')


#' # Get all stats

#' ## Per tool val rates => stats
#' ### RT-qPCR precision
perc_val %>%
  filter(count_group == 'count ≥ 5') %>%
  pull(perc_qPCR_val) %>%
  quantile()

perc_val %>%
  filter(count_group == 'count < 5') %>%
  pull(perc_qPCR_val) %>%
  quantile()

#' ### RNase R precision
perc_val %>%
  filter(count_group == 'count ≥ 5') %>%
  pull(perc_RR_val) %>%
  quantile()

perc_val %>%
  filter(count_group == 'count < 5') %>%
  pull(perc_RR_val) %>%
  quantile()

#' ### amplicon sequencing precision
perc_val %>%
  filter(count_group == 'count ≥ 5') %>%
  pull(perc_amp_val) %>%
  quantile()

perc_val %>%
  filter(count_group == 'count < 5') %>%
  pull(perc_amp_val) %>%
  quantile()

#' ### compound precision
perc_val %>%
  filter(count_group == 'count ≥ 5') %>%
  pull(perc_compound_val) %>%
  quantile()

perc_val %>%
  filter(count_group == 'count < 5') %>%
  pull(perc_compound_val) %>%
  quantile()

#' ### sensitivity

sens %>%
  filter(!count_group_median == 'count ≥ 5') %>%
  filter(!tool == 'circRNA_finder',
         !tool == 'segemehl') %>%
  pull(sensitivity) %>%
  quantile()

sens %>%
  filter(!count_group_median == 'count < 5') %>%
  pull(sensitivity) %>%
  quantile()


#' ## Overall validation rates (ignoring tools) - all circ
#' ### RT-qPCR
cq %>%
  select(circ_id, cell_line, qPCR_val) %>% unique() %>%
  count(qPCR_val)

round(1468/(1468+37), digits = 3)

#' ### RNase R
cq %>%
  select(circ_id, cell_line, RR_val) %>% unique() %>%
  count(RR_val)

round(1308/(1308+85), digits = 3)
round(112/(1308+85+112), digits = 3)

#' ### amplicon sequencing
cq %>%
  select(circ_id, cell_line, amp_seq_val) %>% unique() %>%
  count(amp_seq_val)

cq %>%
  select(circ_id, cell_line, amp_seq_val_detail) %>% unique() %>%
  count(amp_seq_val_detail)

round(1006/(1006+141), digits = 3)
round(335/(1006+141+335), digits = 3)

#' ## Overall validation rates (ignoring tools) - low-abundant circ
#' ### RT-qPCR
cq %>%
  filter(count_group == 'count < 5') %>%
  select(circ_id, cell_line, qPCR_val) %>% unique() %>%
  count(qPCR_val)

round(243/(243+17), digits = 3)

#' ### RNase R
cq %>%
  filter(count_group == 'count < 5') %>%
  select(circ_id, cell_line, RR_val) %>% unique() %>%
  count(RR_val)

round(165/(165+25), digits = 3)

#' ### amplicon sequencing
cq %>%
  filter(count_group == 'count < 5') %>%
  select(circ_id, cell_line, amp_seq_val) %>% unique() %>%
  count(amp_seq_val)

round((63+136)/(63+136+61), digits = 3)

round(136/(63+136), digits = 3)


#' ## Calculate perc of circRNAs upsetplot

#' nr of unique circ that pass all val methods
cq %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line) %>%
  unique() %>%
  filter(qPCR_val == 'pass', RR_val == 'pass', amp_seq_val == 'pass')

round(949/1505, digits = 3)

#' nr of unique circ that fail all val methods
cq %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line) %>%
  unique() %>%
  filter(qPCR_val == 'fail', RR_val == 'fail', amp_seq_val == 'fail')

round(18/1505, digits = 3)


#' nr of unique circ that fail 1 or 2 val methods

cq %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line) %>%
  filter(!is.na(RR_val), !is.na(amp_seq_val)) %>%
  unique() %>%
  mutate(val_nr = str_count(paste(qPCR_val, RR_val, amp_seq_val),
                              pattern = 'pass')) %>%
  #filter(val_nr == 0)
  #filter(val_nr == 3)
  filter(val_nr < 3, val_nr > 0)

round(127/1505, digits = 3)


#' nr of cicrc that have at least one NA
cq %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line) %>%
  unique() %>%
  mutate(NA_nr = str_count(paste(qPCR_val, RR_val, amp_seq_val),
                            pattern = 'NA')) %>%
  #filter(NA_nr == 0)
  #filter(NA_nr == 3)
  filter(NA_nr > 0)

round(411/1505, digits = 3)


#' nr of unique circ that have no NAs
cq %>% 
  filter(!is.na(RR_val), !is.na(amp_seq_val)) %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, cell_line)  %>%
  unique()

#' => recalculate previous numbers based on no NAs
#' all pass
round(949/1094, digits = 3)
#' all fail
round(18/1094, digits = 3)
#' one or two fail
round(127/1094, digits = 3)


