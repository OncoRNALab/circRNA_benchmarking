#' ---
#' title: "Calculations for paper"
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

#' ## Load standard libraries and resolve conflicts
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')


#' ## Read data
all_circ = read_tsv("../data/Supplementary_Table_2_all_circRNAs.txt")

all_circ

cq = read_tsv('../data/Supplementary_Table_3_selected_circRNAs.txt')

sens = read_tsv('../data/Supplementary_Table_6B_sensitivity_values.txt')


#' # How many transcripts have circRNAs?
#' ENST
all_circ %>%
  filter(!is.na(ENST), !ENST == 'ambiguous') %>%
  select(ENST) %>% unique() %>%
  nrow()

# NA or ambiguous => strand needs to be taken into account
all_circ %>% select(circ_id_strand, ENST) %>% unique() %>%
  filter(is.na(ENST) | ENST == 'ambiguous') %>%
  count(ENST)

# total nr

ensembl = read_tsv('/Users/mivromma/Documents/PhD/indexes/Homo_sapiens.GRCh38.103.gtf',
                   col_names = c('chr', 'havana', 'type', 'start', 'end', 'dot', 'strand', 'dot2', 'info'))

ensembl %>% filter(type == "transcript") %>%
  separate(info, into = c('gene_id', 'gene_version', 'transcrip_id', 'rest'), 
           extra = 'merge', sep = ';') %>%
  select(transcrip_id) %>% unique()

17461/234393

can = read_tsv('known_exons_GRCh38.103_canonical.bed', col_names = c('chr', 'start', 'end', 'exon', 'dot', 'strand'))

can %>% separate(exon, into = c('ENST', 'ENST_nr', 'exon', 'exon_nr'), sep = '_') %>%
  select(ENST) %>% unique() %>%
  inner_join(all_circ %>% select(ENST) %>% unique())

17461/60427

#' # CircRNA length stats
#' ### All circRNAs
#' with introns
all_circ %>% select(chr, start, end, estim_len_in) %>% unique() %>%
  pull(estim_len_in) %>% quantile()

#' no introns
all_circ %>% select(chr, start, end, estim_len_ex) %>% unique() %>%
  filter(!(is.na(estim_len_ex) | estim_len_ex == 'ambiguous')) %>%
  mutate(estim_len_ex = as.numeric(estim_len_ex)) %>%
  pull(estim_len_ex) %>% quantile()

#' ### CircRNAs that don't get validated with RNase R

#' do not get validated (no introns)
cq %>% filter(qPCR_val == 'pass',
              RR_val == 'fail',
              amp_seq_val == 'pass') %>%
  select(chr, start, end, estim_len_ex) %>% unique() %>%
  filter(!(is.na(estim_len_ex) | estim_len_ex == 'ambiguous')) %>%
  mutate(estim_len_ex = as.numeric(estim_len_ex)) %>%
  pull(estim_len_ex) %>% quantile()


#' # BSJ count distribution

all_circ %>% 
  filter(!tool == "Sailfish-cir") %>%
  pull(BSJ_count) %>%
  quantile()

#' how many with BSJ < 5
all_circ %>% 
  filter(!tool == "Sailfish-cir") %>%
  count(count_group)

round(945201/(945201+146700), 3)

#' how many with BSJ >= 2
all_circ %>% 
  filter(!tool == "Sailfish-cir") %>%
  filter(BSJ_count >= 2) %>% nrow()

round(503237/(945201+146700), 3)

#' # How many circRNAs are detected on different strands
all_circ %>% 
  filter(!strand == 'unknown') %>%
  select(chr, start, end) %>% unique() %>%
  nrow()

all_circ %>% 
  filter(!strand == 'unknown') %>%
  select(chr, start, end, strand) %>% unique() %>%
  nrow()

round((190314-174009)/174009, 3)


#' # CircRNA numbers

#' nr of circ detected per tool
all_circ %>%
  filter(cell_line == 'HLF') %>%
  group_by(tool) %>%
  count() %>%
  arrange(n)

#' nr of unique circ
all_circ %>% select(circ_id) %>%
  unique()

#' total nr of circ
nrow(all_circ)
