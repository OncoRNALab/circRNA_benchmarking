
#' ---
#' title: "Fix bug in data due to position shift"
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


#' # read in data

all_circ = read_tsv('../data/Supplementary_Table_2_all_circRNAs.txt')

all_circ

corrected_circ = read_tsv('02_corrected_circ.txt') %>%
  select(circ_id, pos_shift) %>%
  unique() %>%
  rename(new_circ_id = circ_id,
         circ_id = pos_shift)

corrected_circ

#' # fix selected circ table

cq = read_tsv('data_with_mistake/selected_circ_cq.txt')

cq %>% select(circ_id, cell_line) %>% unique()

#' ## first fix those spec circ

cq_fix = cq %>%
  inner_join(corrected_circ) %>%
  select(tool, cell_line, FWD_primer, REV_primer, FWD_pos, FWD_length, REV_pos, REV_length,
         FWD_TM, REV_TM, FWD_GC, REV_GC, amplicon, primer_unique, BSJ_count, Cq_min_untreated, Cq_max_untreated,
         Cq_min_treated, Cq_max_treated, new_circ_id) %>%
  rename(circ_id = new_circ_id) 

cq_fix

# add anno and fix primer positions

cq_fix = cq_fix %>%
  left_join(all_circ %>% 
              select(chr, start, end, strand, BSJ_count, cell_line, tool, circ_id,
                     circ_id_strand, count_group)) %>%
  mutate(FWD_pos = FWD_pos - 1,
         REV_pos = REV_pos + 1)
cq_fix

#' add to original dataframe

cq = cq %>%
  anti_join(corrected_circ) %>%
  bind_rows(cq_fix)

cq

cq %>% select(circ_id, cell_line) %>% unique()


#' ## fix other mistake from repeated circ
#' there are 11 circ that are accidentally measured twice (same primers) => keep the lowest Cq values
cq %>%
  select(circ_id, cell_line, Cq_min_untreated, Cq_max_untreated,
         Cq_min_treated, Cq_max_treated, FWD_primer, REV_primer) %>%
  unique() %>%
  group_by(circ_id, cell_line) %>% 
  filter(n() > 1)

cq_fix = cq %>%
  select(circ_id, cell_line, Cq_min_untreated, Cq_max_untreated,
         Cq_min_treated, Cq_max_treated,FWD_primer, 
         REV_primer, FWD_pos, FWD_length, REV_pos, REV_length, FWD_TM, REV_TM,
         FWD_GC, REV_GC, amplicon, primer_unique) %>%
  unique() %>%
  group_by(circ_id, cell_line) %>% 
  filter(n() > 1) %>%
  mutate(tmp_cq = (Cq_min_untreated + Cq_max_untreated)/2) %>%
  slice(which.min(tmp_cq)) %>%
  ungroup() %>%
  select(-tmp_cq)


cq_fix

# add back tool annotation
cq_fix = cq_fix %>%
  inner_join(cq %>%
               select(chr, start, end, strand, circ_id, circ_id_strand, tool, cell_line,
                      BSJ_count, count_group))

cq = cq %>%
  anti_join(cq_fix %>% select(circ_id, cell_line) %>% unique()) %>%
  bind_rows(cq_fix)


#' ## fix primer unique label
cq %>% count(primer_unique)


cq_fix = cq %>%
  select(circ_id, FWD_primer, REV_primer) %>%
  unique() %>%
  group_by(FWD_primer, REV_primer) %>%
  mutate(primer_unique = ifelse(n() > 1, 'no', 'yes')) %>%
  ungroup() %>% unique()

cq_fix

cq = cq %>%
  select(-primer_unique) %>%
  left_join(cq_fix)

cq %>% count(primer_unique)


#' put in right order again and save data_frame

cq = cq %>%
  select(chr, start, end, strand, circ_id, circ_id_strand, tool, cell_line, FWD_primer, 
         REV_primer, FWD_pos, FWD_length, REV_pos, REV_length, FWD_TM, REV_TM,
         FWD_GC, REV_GC, amplicon, primer_unique, BSJ_count, count_group, Cq_min_untreated,
         Cq_max_untreated, Cq_min_treated, Cq_max_treated)

cq %>% write_tsv('../data/details/selected_circ_cq.txt')


#' #  fix amp seq file

amp = read_tsv('data_with_mistake/amp_seq.txt')

amp

#' ## first fix specific circ

amp_fix = amp %>% 
  inner_join(corrected_circ) %>%
  mutate(circ_id = new_circ_id) %>%
  select(-new_circ_id)

amp_fix

amp = amp %>%
  anti_join(corrected_circ) %>%
  bind_rows(amp_fix)

amp

#' ## then fix doubles

amp %>% 
  unique() %>% # there are some exact doubles, but this is not a big problem
  group_by(circ_id, cell_line) %>%
  filter(n() > 1)

amp = amp %>%
  filter(!(circ_id == 'chr1:147252622-147259918' & amp_seq_val == "fail")) %>%
  unique()

#' write new dataframe

amp %>% write_tsv('../data/details/amp_seq.txt')


#' # fix Sup T4

table4 = read_tsv('../data/Supplementary_Table_4_all_circRNAs_treated.txt')

table4 %>% inner_join(corrected_circ) %>% view()

table4 %>% filter(strand == '+',
                  tool %in% c('KNIFE', 'NCLscan', 'NCLcomparator'))
