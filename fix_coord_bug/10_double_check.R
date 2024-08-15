#' ---
#' title: "Double chek all tables"
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


#' Sup table 2

sup = read_tsv('../data/Supplementary_Table_2_all_circRNAs.txt')
sup = read_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')
sup = read_tsv('../data/Supplementary_Table_4_all_circRNAs_treated.txt')
sup = read_tsv('../data/Supplementary_Table_5_RNase_R_enrichment_seq.txt.gz')


sup
  
# should be: chr 18 8718423 8720496
sup %>% filter(circ_id == 'chr18:8718423-8720496') %>% pull(tool) %>% unique()
# but is in final frame: 8718424 8720495
sup %>% filter(circ_id == 'chr18:8718424-8720495') %>% pull(tool) %>% unique()



sup %>% 
  filter(circ_id == 'chr14:92126121-92143294') %>% view()
sup %>% 
  filter(circ_id == 'chr1:147252622-147259918') %>% view()

sup %>% 
  filter(circ_id == 'chr3:93995873-94003908') %>% view()


sup %>% 
  filter(circ_id == 'chr12:27714779-27724186') %>% view()
sup %>% 
  filter(circ_id == 'chrX:47243395-47245482') %>% view()
sup %>% 
  filter(circ_id == 'chr5:14270824-14287063') %>% view()
sup %>% filter(is.na(n_detected))

sup %>% select(circ_id) %>% unique()

sup = read_tsv('../data/Supplementary_Table_2_all_circRNAs.txt')

sup %>% select(circ_id_strand, ENST) %>% unique() %>%
  filter(!is.na(ENST),
         !ENST == 'ambiguous') %>%
  count(ENST)

201635/385378
39326/385378
144417/385378
(39326+144417)/385378

sup %>%
  filter(!is.na(ENST), !ENST == 'ambiguous') %>%
  select(ENST) %>% unique() %>%
  nrow()

sup %>% select(circ_id) %>% unique()

sup %>% select(circ_id, n_detected) %>% unique() %>%
  count(n_detected == 1)

193759/(193759+161011)
180912/(180912+152281)

sup %>% select(circ_id, n_detected, cell_line) %>% unique() %>%
  filter(n_detected == 16) %>%
  select(circ_id) %>% unique()

sup = read_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')

sup %>% filter(compound_val == 'pass') %>% select(circ_id) %>% unique()


sup %>% select(circ_id, cell_line) %>% unique() %>%
  group_by(circ_id) %>% summarise(n_cell = n()) %>% filter(n_cell > 1) %>%
  count(n_cell)

sup %>% filter(!count_group == 'count < 5',
               !is.na(amp_seq_val)) %>% count(tool) %>% pull(n) %>% quantile
sup %>% count(RR_val_detail)

sup %>% select(circ_id, qPCR_val, RR_val, amp_seq_val, amp_seq_val_detail) %>% unique() %>%
  filter(qPCR_val == 'pass',
         RR_val == 'fail',
         amp_seq_val == 'fail')

sup %>%
  filter(n_detected > 1) %>%
  select(circ_id, qPCR_val, RR_val, amp_seq_val, compound_val) %>%
  unique() %>%
  filter(qPCR_val == 'fail' | RR_val == 'fail' | amp_seq_val == 'fail')

sup %>% select(circ_id) %>% unique()

sup = read_tsv('../data/Supplementary_Table_6_tool_ranking.txt')

sup %>% filter(count_group == 'count ≥ 5') %>%
  filter(compound_precision > 0.9)


# fix tool annotation

t9 = read_tsv('../data/Supplementary_Table_9_top_tool_combinations.txt')

t9 %>%
  count(splicing_1)
t9 %>%
  count(splicing_2)

t9 = t9 %>%
  mutate(splicing_2 = ifelse(tool_2 == 'NCLcomparator',
                             'canonical',
                             splicing_2))
t9 %>% write_tsv('../data/Supplementary_Table_9_top_tool_combinations.txt')


sup2 %>% filter(circ_id == 'chr1:90937484-90982370') %>% select(tool) %>% unique()


# check numbers figure 21
val = read_tsv('../data/Supplementary_Table_6A_precision_values.txt')
val %>%
  filter(count_group == 'count ≥ 5') %>%
  mutate(per_pass = nr_amp_val/80,
         per_fail = nr_amp_fail/80) %>% view()
