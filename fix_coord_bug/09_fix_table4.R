
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

table4 = read_tsv('data_with_mistake/Supplementary_Table_4_all_circRNAs_treated.txt')

#' # fix Sup T4
#' follow the same logical as script 01_fix_coord_bug.R and 03_add_annotation.R (but no detailed transcript annotation for now)


table4 %>% filter(strand == '+',
                  tool %in% c('KNIFE', 'NCLscan', 'NCLcomparator'))

# should be: chr 18 8718423 8720496
table4 %>% filter(circ_id == 'chr18:8718423-8720496') %>% pull(tool) %>% unique()
# but is in final frame: 8718424 8720495
table4 %>% filter(circ_id == 'chr18:8718424-8720495') %>% pull(tool) %>% unique()


#' fix positions
corrected_circ = table4 %>%
  filter(tool %in% c("KNIFE","NCLscan","NCLcomparator"),
         strand == '+') %>%  # this selects all the circRNAs that were wrongly 'corrected'
  mutate(pos_shift = circ_id, # to keep original circ_id
         end = end + 1,  # because these have been wrongfully been -1
         start = start -1) %>% # because this is the actual conversion that was necessary
  mutate(circ_id = paste(chr, ':', start, "-", end, sep = ""),
         circ_id_strand = paste(chr, ':', start, "-", end, "/", strand, sep = "")) %>%
  select(chr, start, end, strand, BSJ_count, cell_line, RNaseR, tool, circ_id, 
         circ_id_strand, count_group, pos_shift) # only select columns that are not problematic

corrected_circ


# should be: chr 18 8718423 8720496
corrected_circ %>% filter(circ_id == 'chr18:8718423-8720496') %>% pull(tool) %>% unique()
# but is in final frame: 8718424 8720495
corrected_circ %>% filter(circ_id == 'chr18:8718424-8720495') %>% pull(tool) %>% unique()

#' add all annotation to non-unique circRNAs

corrected_circ =
  corrected_circ %>%
  left_join(table4 %>% 
              select(circ_id_strand, estim_len_in, ENST, estim_len_ex,
                     nr_exons, start_match, end_match, ss_motif) %>%
              unique())




#' put in the same order
table4_new = table4 %>%
  filter(!(strand == '+' &
           tool %in% c('KNIFE', 'NCLscan', 'NCLcomparator'))) %>%
  bind_rows(corrected_circ) %>%
  select(chr, start, end, strand, BSJ_count, cell_line, RNaseR, tool, circ_id, circ_id_strand, 
         count_group, estim_len_in, ss_motif, estim_len_ex, ENST,
         nr_exons, start_match, end_match)


table4_new

# should be: chr 18 8718423 8720496
table4_new %>% filter(circ_id == 'chr18:8718423-8720496') %>% pull(tool) %>% unique()
# but is in final frame: 8718424 8720495
table4_new %>% filter(circ_id == 'chr18:8718424-8720495') %>% pull(tool) %>% unique()

#' write new file

table4_new %>%
  write_tsv("../data/Supplementary_Table_4_all_circRNAs_treated.txt")
