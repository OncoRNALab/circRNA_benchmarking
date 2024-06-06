
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

#' # analyse bug
#' these tools were correctly assigned as 1-based tools: "circtools","CIRI2","CIRIquant","KNIFE","NCLscan","NCLscan_post","sailfish-cir"
#' to fix this, the start position was replaced with 'start -1'
#' later, it became clear that the circRNAs on the + strand of "KNIFE","NCLscan", and "NCLscan_post" are reported as 'chr:end-start' which lead to the wrong 'correction' from 1-based to 0-based

# example
# original entry: chr18|MTCL1:8720496|MTCL1:8718424|rev|+	894	1	0	RNA015434_S1 (1-based, start > end)
# chr18   8720495   8718424 +        894 HLF       untreated KNIFE chr18:8720495-8718424    
# should be: chr 18 8718423 8720496
all_circ %>% filter(circ_id == 'chr18:8718423-8720496') %>% pull(tool) %>% unique()
# but is in final frame: 8718424 8720495
all_circ %>% filter(circ_id == 'chr18:8718424-8720495') %>% pull(tool) %>% unique()


#' # fix bug

corrected_circ = all_circ %>%
  filter(tool %in% c("KNIFE","NCLscan","NCLcomparator"),
         strand == '+') %>%  # this selects all the circRNAs that were wrongly 'corrected'
  filter(end > start) %>%
  mutate(pos_shift = circ_id, # to keep original circ_id
         end = end + 1,  # because these have been wrongfully been -1
         start = start -1) %>% # because this is the actual conversion that was necessary
  mutate(circ_id = paste(chr, ':', start, "-", end, sep = ""),
         circ_id_strand = paste(chr, ':', start, "-", end, "/", strand, sep = "")) %>%
  select(chr, start, end, strand, BSJ_count, cell_line, tool, circ_id, 
         circ_id_strand, count_group, pos_shift) # only select columns that are not problematic

corrected_circ

corrected_circ %>% write_tsv('02_corrected_circ.txt')


# should be: chr 18 8718423 8720496
corrected_circ %>% filter(circ_id == 'chr18:8718423-8720496') %>% pull(tool) %>% unique()
# but is in final frame: 8718424 8720495
corrected_circ %>% filter(circ_id == 'chr18:8718424-8720495') %>% pull(tool) %>% unique()


#' # check effect on validation

cq = read_tsv('../data/details/selected_circ_cq.txt')
cq

cq_affected = cq %>% 
  inner_join(corrected_circ %>% 
               select(pos_shift, circ_id) %>%
               mutate(new_id = circ_id,
                      circ_id = pos_shift) %>%
               unique())

cq_affected             

cq_affected %>%
  count(tool)




