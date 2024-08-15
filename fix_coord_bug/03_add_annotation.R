
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

all_circ = read_tsv('../data/Supplementary_Table_2_all_circRNAs.txt')#' # add known annotation
corrected_circ = read_tsv('02_corrected_circ.txt')


#' add info to non-unique circRNAs

corrected_circ =
  corrected_circ %>%
  left_join(all_circ %>% 
              select(circ_id_strand, n_db, dbs, estim_len_in, ENST, estim_len_ex,
                     nr_exons, start_match, end_match, ss_motif) %>%
              unique())

#' make subset of those that need annotation
sub_circ = corrected_circ %>% 
  anti_join(all_circ %>% 
              select(circ_id_strand) %>%
              unique()) %>% 
  select(chr, start, end, strand, circ_id, circ_id_strand,
         pos_shift) %>%
  unique()

sub_circ

#' ## add database information

#' # comparison databases

db_circ = read_tsv('../data/details/circ_db_hg38.txt')

db_circ

db_circ = db_circ %>%
  group_by(circ_id) %>%
  summarize(n_db = n(),
            dbs = paste(database, collapse = '/'))

db_circ

sub_circ = sub_circ %>%
  left_join(db_circ)

sub_circ

#' check if similar to all circ

sub_circ %>% count(n_db)
all_circ %>% count(n_db)

#' # add transcript info
#' ## estimated length with introns

sub_circ = sub_circ %>% 
  mutate(estim_len_in = end - start)

sub_circ

#' ## match circ info with transcript info

#' make bed file
sub_circ %>% 
  mutate(score = 0) %>%
  select(chr, start, end, circ_id_strand, score, strand) %>% 
  write_tsv('04_sub_circ_in.bed', col_names = F)

#' commands for bed tools
#' docker run -v ${PWD}:/usr/data -it --entrypoint /bin/bash quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0
#' cd /usr/data
#' bedtools sort -i 04_sub_circ_in.bed > 05_sub_circ_in.sorted.bed
#' bedtools intersect -wa -wb -s -a 05_sub_circ_in.sorted.bed -b known_exons_GRCh38.103_canonical.bed > 06_sub_circ_canonical_match.bed

#' read in BED intersect results
exons_long = read_tsv('06_sub_circ_canonical_match.bed',
                  col_names = c('overlap_chr', 'overlap_start', 'overlap_end', 'circ_id', 'score',
                                'overlap_strand', 'exon_chr', 'exon_start', 'exon_end', 'exon', 'exon_score',
                                'exon_strand', 'overlap'))
exons_long

#' clean up data frame
exons_long = exons_long %>%
  separate(circ_id, into = c('chr', 'pos'), remove = F, sep = ":") %>%
  separate(pos, into = c('start', 'end'), sep = '-', extra = "merge") %>%
  separate(end, into = c('end', 'strand'), sep = "/") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end)) %>%
  separate(exon, into = c('ENST', 'ENST_v', 'skip', 'exon_nr'), sep = '_', remove = F)

exons_long

#' get unique ENST for each circ or ambiguous + estimated length without introns

exons = exons_long %>%
  group_by(circ_id) %>%
  filter(length(unique(ENST)) == 1) %>%
  summarize(ENST = min(ENST), 
            estim_len_ex = sum(overlap),
            nr_exons = n()) %>%
  separate(circ_id, into = c('circ_id', 'strand'), sep = '/') %>% 
  mutate(estim_len_ex = as.character(estim_len_ex),
         nr_exons = as.character(nr_exons)) %>%
  bind_rows(exons %>%
              group_by(circ_id) %>%
              filter(length(unique(ENST)) > 1) %>%
              select(circ_id) %>%
              separate(circ_id, into = c('circ_id', 'strand'), sep = '/') %>%
              mutate(ENST = "ambiguous", estim_len_ex = "ambiguous", nr_exons = "ambiguous") %>%
              unique())

exons

#' add info to dataframe and add NA if no match
sub_circ = 
  sub_circ %>%
  left_join(exons)

sub_circ

#' ## exon match info
# make df with start and end info

exons_se = exons %>% 
  left_join(exons_long %>% 
              select(circ_id, ENST, start, exon_start, end, exon_end, exon) %>%
              separate(circ_id, into = c('circ_id', 'strand'), sep = '/'))

exons_se

BSJ_match = exons_se %>% 
  filter(!is.na(exon_start & !is.na(exon_end))) %>%
  mutate(start_match = ifelse(start == exon_start, exon, "no_match")) %>%
  filter(!start_match == 'no_match') %>%
  select(circ_id, strand, ENST, estim_len_ex, nr_exons, start, end, start_match) %>%
  full_join(exons_se %>% 
              filter(!is.na(exon_start & !is.na(exon_end))) %>%
              mutate(end_match = ifelse(end == exon_end, exon, "no_match")) %>%
              filter(!end_match == 'no_match') %>%
              select(circ_id, strand, ENST, estim_len_ex, nr_exons, start, end, end_match))

BSJ_match

#' add to dataframe

sub_circ = 
  sub_circ %>%
  left_join(BSJ_match)

sub_circ


#' # get splicing sequence with fastahack

#' make files with SS positions
SS_file = sub_circ %>%
  select(chr, start, end) %>%
  unique()

SS_file

#' for SS acceptor (2 nt before circ start position)
SS_file %>%
  # only -1 because of simultanous 0-based to 1-based annotation
  select(chr, start) %>% mutate(end = start, start = start - 1) %>% unique() %>%
  mutate(ss = paste(str_remove(chr, 'chr'), ":",start,"-", end, sep = '')) %>%
  select(ss) %>%
  write_tsv("07_circ_in_acc.txt", col_names = F)

#' for SS donor (2 nt after circ end position)
SS_file %>%
  # 0-based to 1-based annotation
  select(chr, end) %>% mutate(start = end + 1, end = end + 2) %>% unique() %>%
  mutate(ss = paste(str_remove(chr, 'chr'), ":",start,"-", end, sep = '')) %>%
  select(ss) %>%
  write_tsv("07_circ_in_donor.txt", col_names = F)


#' fastahack commands
#' docker run -v ${PWD}:/usr/data -it --entrypoint oncornalab/primerxl_circ:v0.30
#' cd /usr/data
#' /bin/fastahack-1.0.0/fastahack -i Homo_sapiens.GRCh38.dna.primary_assembly.fa
#' cat 07_circ_in_acc.txt | /bin/fastahack-1.0.0/fastahack -c Homo_sapiens.GRCh38.dna.primary_assembly.fa > 08_circ_out_acc.txt
#' cat 07_circ_in_donor.txt | /bin/fastahack-1.0.0/fastahack -c Homo_sapiens.GRCh38.dna.primary_assembly.fa > 08_circ_out_donor.txt
#' 
#' first all for acceptor
circ_ss_acc = 
  read_tsv("07_circ_in_acc.txt", col_names = c("ss_id")) %>%
  bind_cols(read_tsv("08_circ_out_acc.txt", col_names = c('ss_seq'))) %>%
  separate(ss_id, into = c("chr", "tmp"), sep = ':') %>%
  separate(tmp, into = c("start", 'end'), sep= "-") %>%
  mutate(start = as.numeric(end), chr = paste("chr", chr, sep = "")) %>% select(-end) # SS info to circRNA info

circ_ss_acc

# then all for donor
circ_ss_donor = 
  read_tsv("07_circ_in_donor.txt", col_names = c("ss_id")) %>%
  bind_cols(read_tsv("08_circ_out_donor.txt", col_names = c('ss_seq'))) %>%
  separate(ss_id, into = c("chr", "tmp"), sep = ':') %>%
  separate(tmp, into = c("start", 'end'), sep= "-") %>%
  mutate(end = as.numeric(start) - 1, chr = paste("chr", chr, sep = "")) %>% select(-start) # SS info to circRNA info


circ_ss_donor

#' link fastahack info back to original circRNAs

sub_circ = sub_circ %>% 
  left_join(circ_ss_donor %>% rename(ss_donor = ss_seq)) %>%
  left_join(circ_ss_acc %>% rename(ss_acceptor = ss_seq))

sub_circ

#' fix for neg strand not needed because all are on pos strand

#' paste acceptor and donor together
sub_circ = sub_circ %>%
  mutate(ss_motif = paste(ss_acceptor, "N", ss_donor, sep = ""))


#' # put new info with all_circ
#' add sub_circ and clean up corrected circ dataframe

corrected_circ_merged = 
  # all the ones that already had annotation (non-unique)
  corrected_circ %>%
  inner_join(all_circ %>% 
              select(circ_id_strand) %>%
              unique()) %>%
  select(-pos_shift) %>%
  # all the ones that needed new annotation (non-unique)
  bind_rows(corrected_circ %>% 
              anti_join(all_circ %>% 
                          select(circ_id_strand) %>%
                          unique()) %>% 
              select(chr, start, end, strand, BSJ_count, cell_line, tool, circ_id, circ_id_strand,
                     count_group) %>%
              # add new annotation
              left_join(sub_circ %>%
                          select(-pos_shift, -ss_donor, -ss_acceptor)))

corrected_circ_merged

#' merge with all circ dataframe
all_circ_new = all_circ %>% 
  filter(!(tool %in% c("KNIFE","NCLscan","NCLcomparator") &
             strand == '+')) %>%
  bind_rows(corrected_circ_merged)


#' double check
# should be: chr 18 8718423 8720496
all_circ %>% filter(circ_id == 'chr18:8718423-8720496') %>% pull(tool) %>% unique()
all_circ_new %>% filter(circ_id == 'chr18:8718423-8720496') %>% pull(tool) %>% unique()

# but is in final frame: 8718424 8720495
all_circ %>% filter(circ_id == 'chr18:8718424-8720495') %>% pull(tool) %>% unique()
all_circ_new %>% filter(circ_id == 'chr18:8718424-8720495') %>% pull(tool) %>% unique()


#' # update nr of tools

n_tools =
  all_circ_new %>%
  select(circ_id, tool, cell_line)%>% unique() %>%
  group_by(circ_id, cell_line) %>%
  summarise(
    n_detected = n(),
    tools = paste(tool, collapse = '/'))

n_tools

all_circ_new =
  all_circ_new %>%
  select(-n_detected, -tools) %>%
  left_join(n_tools)

all_circ_new %>% select(circ_id, n_detected) %>% unique() %>% filter(n_detected == 16)
all_circ_new %>% filter(n_detected == 16) %>% select(circ_id) %>% unique()

#' # update BSJ count median + group

all_circ_new = all_circ_new %>%
  group_by(circ_id, cell_line) %>%
  mutate(BSJ_count_median = median(BSJ_count),
         count_group_median = ifelse(BSJ_count_median >= 5, 'count ≥ 5', 'count < 5')) %>%
  ungroup()
all_circ_new

#' put in the same order
all_circ_new = all_circ_new %>%
  select(chr, start, end, strand, BSJ_count, cell_line, tool, circ_id, circ_id_strand, 
         count_group, n_detected, tools, n_db, dbs, estim_len_in, ENST, estim_len_ex,
         nr_exons, start_match, end_match, ss_motif, BSJ_count_median, count_group_median)


#' write new file

all_circ_new %>%
  write_tsv("../data/Supplementary_Table_2_all_circRNAs.txt")
  