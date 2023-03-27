#' ---
#' title: "Combination of tools"
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

#' ## Loads standard libraries and resolve conflicts
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')

#' ## Load specific libraries
library(gplots)
library(ggpubr)
library("quantreg") 
library(ggrepel)


#' ## Read data
all_circ = read_tsv('../data/Supplementary_Table_2_all_circRNAs.txt')
val = read_tsv('../data/Supplementary_Table_6A_precision_values.txt')

all_circ
val


#' # Two tools
#' only for circRNAs in BSJ count group ≥ 5
combo = tibble()


for (tool_1 in all_circ %>% pull(tool) %>% unique()) {
  for (tool_2 in all_circ %>% pull(tool) %>% unique()) {
    for (cell in all_circ %>% pull(cell_line) %>% unique()) {
      nr_union = all_circ %>%
        filter(tool %in% c(tool_1, tool_2),
               count_group == 'count ≥ 5',
               cell_line == cell) %>%
        select(circ_id) %>%
        unique() %>% nrow()

      nr_intersection = all_circ %>%
        filter(tool == tool_1,
               count_group == 'count ≥ 5',
               cell_line == cell) %>%
        select(circ_id) %>%
        inner_join(all_circ %>%
                     filter(tool == tool_2,
                            count_group == 'count ≥ 5',
                            cell_line == cell) %>%
                     select(circ_id)) %>%
        unique() %>% nrow()

      combo = combo %>%
        bind_rows(tibble(tool_1, tool_2, cell, nr_union, nr_intersection))

    }
  }
}

combo

#' clean up dataframes and add validation rates
combo = combo %>%
  rename(cell_line = cell) %>%
  left_join(val %>%
              filter(count_group == 'count ≥ 5') %>%
              rename(tool_1 = tool, perc_compound_val_1 = perc_compound_val) %>%
              select(tool_1, perc_compound_val_1, count_group)) %>%
  left_join(val %>%
              filter(count_group == 'count ≥ 5') %>%
              rename(tool_2 = tool, perc_compound_val_2 = perc_compound_val) %>%
              select(tool_2, perc_compound_val_2, count_group))

#' add total nr of circ per tool
combo = combo %>%
  left_join(all_circ %>%
              filter(count_group == 'count ≥ 5') %>%
              group_by(cell_line, tool) %>%
              select(circ_id, cell_line, tool) %>%
              unique() %>%
              summarise(total_n = n()) %>% ungroup() %>%
              rename(tool_1 = tool, total_n_1 = total_n)) %>%
  left_join(all_circ %>%
              filter(count_group == 'count ≥ 5') %>%
              group_by(cell_line, tool) %>%
              select(circ_id, cell_line, tool) %>%
              unique() %>%
              summarise(total_n = n()) %>% ungroup() %>%
              rename(tool_2 = tool, total_n_2 = total_n)) %>%
  mutate(w_val_rate = ((perc_compound_val_1 * total_n_1) + (perc_compound_val_2 * total_n_2)) /
           (total_n_1 + total_n_2)) %>%
  filter(!is.na(perc_compound_val_1), !is.na(perc_compound_val_2))

combo

#' add tool letters
combo = combo %>%
  left_join(read_tsv('../data/details/tool_annotation.txt') %>%
              select(tool, tool_lt) %>%
              rename(tool_1 = tool, tool_lt_1 = tool_lt)) %>%
  left_join(read_tsv('../data/details/tool_annotation.txt') %>%
              select(tool, tool_lt) %>%
              rename(tool_2 = tool, tool_lt_2 = tool_lt)) %>%
  mutate(tool_combo = paste(tool_lt_1, tool_lt_2, sep = ""),
         tool_combo = map_chr(str_split(tool_combo, ""), ~str_c(str_sort(unique(.x)), collapse = "")))

combo


#' add total nr of circ for that cell line and calculate perc increase
total_cl = all_circ %>%
  group_by(cell_line) %>%
  filter(count_group == 'count ≥ 5') %>%
  select(circ_id, cell_line) %>%
  unique() %>%
  count(cell_line) %>%
  rename(total_cell_line = n) %>% ungroup()

total_cl

combo = combo %>%
  left_join(total_cl) %>%
  mutate(perc_increase_t1  = (nr_union - total_n_1)/total_n_1)

combo

#' save as a dataframe for Supplementary table 5

combo %>% write_tsv('../data/Supplementary_Table_7_combo_2tools.txt')


#' ## Generate a top 10 list of combinations (Supplementary Table 9)

#'  remove the combinations of twice the same tool and only have one line per combo
combo_top = combo %>%
  filter(!tool_1 == tool_2) %>%
  filter(tool_1 < tool_2)

#' take all cell lines together and save as separate columns
combo_top = combo_top %>%
  select(tool_1, tool_2, nr_union, cell_line) %>%
  pivot_wider(names_from = cell_line, values_from = nr_union) %>%
  rename(nr_union_HLF = HLF, nr_union_NCI_H23 = `NCI-H23`, nr_union_SW480 = SW480) %>%
  left_join(combo_top %>%
              select(tool_1, tool_2, w_val_rate, cell_line) %>%
              pivot_wider(names_from = cell_line, values_from = w_val_rate) %>%
              rename(w_precision_HLF = HLF, w_precision_NCI_H23 = `NCI-H23`, w_precision_SW480 = SW480))

#' generate a rank for each metric
combo_top = combo_top %>%
  mutate(w_precision_HLF_rank = dense_rank(desc(w_precision_HLF)),
         nr_union_HLF_rank = dense_rank(desc(nr_union_HLF)),
         w_precision_NCI_H23_rank = dense_rank(desc(w_precision_NCI_H23)),
         nr_union_NCI_H23_rank = dense_rank(desc(nr_union_NCI_H23)),
         w_precision_SW480_rank = dense_rank(desc(w_precision_SW480)),
         nr_union_SW480_rank = dense_rank(desc(nr_union_SW480)))

#' select the top 5 in each category
combo_top = combo_top %>%
  slice_max(w_precision_HLF, n = 5) %>%
  bind_rows(combo_top %>% slice_max(nr_union_HLF, n = 5)) %>%
  bind_rows(combo_top %>% slice_max(w_precision_NCI_H23, n = 5)) %>%
  bind_rows(combo_top %>% slice_max(nr_union_NCI_H23, n = 5)) %>%
  bind_rows(combo_top %>% slice_max(w_precision_SW480, n = 5)) %>%
  bind_rows(combo_top %>% slice_max(nr_union_SW480, n = 5)) %>%
  unique()

combo_top

#' add tool annotation

tool_annotation = read_tsv('/Users/mivromma/Documents/PhD/GitHub/circRNA_benchmarking/data/details/tool_annotation.txt')
tool_annotation

combo_top = combo_top %>%
  left_join(tool_annotation %>%
              rename(tool_1 = tool, approach_1 = approach, lin_annotation_1 = lin_annotation,
                     strand_anno_1 = strand_anno, splicing_1 = splicing, BSJ_filter_1 = BSJ_filter) %>%
              select(-tool_lt)) %>%
  left_join(tool_annotation %>%
              rename(tool_2 = tool, approach_2 = approach, lin_annotation_2 = lin_annotation,
                     strand_anno_2 = strand_anno, splicing_2 = splicing, BSJ_filter_2 = BSJ_filter) %>%
              select(-tool_lt))

combo_top

#' save as sup table 9
combo_top %>%
  select("tool_1", "tool_2", "nr_union_HLF", "nr_union_HLF_rank", "nr_union_NCI_H23", 
         "nr_union_NCI_H23_rank", "nr_union_SW480", "nr_union_SW480_rank" , "w_precision_HLF", 
         "w_precision_HLF_rank",  "w_precision_NCI_H23", "w_precision_NCI_H23_rank",
         "w_precision_SW480", "w_precision_SW480_rank", "approach_1", "approach_2",
         "lin_annotation_1", "lin_annotation_2", "strand_anno_1", "strand_anno_2",
         "splicing_1", "splicing_2", "BSJ_filter_1", 'BSJ_filter_2') %>% 
  arrange(tool_1, tool_2) %>%
  write_tsv('../data/Supplementary_Table_9_top_tool_combinations.txt')


#' # Three tools
#' only for circRNAs in BSJ count group ≥ 5
combo_3 = tibble()

for (cell in all_circ %>% pull(cell_line) %>% unique()) {
  for (tool_1 in all_circ %>% pull(tool) %>% unique()) {
    for (tool_2 in all_circ %>% pull(tool) %>% unique()) {
      for (tool_3 in all_circ %>% pull(tool) %>% unique()) {
        nr_union = all_circ %>%
          filter(count_group == 'count ≥ 5') %>%
          filter(tool %in% c(tool_1, tool_2, tool_3)) %>%
          filter(cell_line == cell) %>%
          select(circ_id) %>%
          unique() %>% nrow()

        combo_3 = combo_3 %>%
          bind_rows(tibble(tool_1, tool_2, tool_3, nr_union, cell))
      }
    }
  }
}

combo_3


#' clean up dataframe and add validation rates
combo_3 = combo_3 %>%
  rename(cell_line = cell) %>%
  left_join(val %>%
              filter(count_group == 'count ≥ 5') %>%
              rename(tool_1 = tool, perc_compound_val_1 = perc_compound_val) %>%
              select(tool_1, perc_compound_val_1, count_group)) %>%
  left_join(val %>%
              filter(count_group == 'count ≥ 5') %>%
              rename(tool_2 = tool, perc_compound_val_2 = perc_compound_val) %>%
              select(tool_2, perc_compound_val_2, count_group)) %>%
  left_join(val %>%
              filter(count_group == 'count ≥ 5') %>%
              rename(tool_3 = tool, perc_compound_val_3 = perc_compound_val) %>%
              select(tool_3, perc_compound_val_3, count_group))

#' add nr of circ from each group

combo_3 = combo_3 %>%
  left_join(all_circ %>%
              filter(count_group == 'count ≥ 5') %>%
              group_by(cell_line, tool) %>%
              select(circ_id, cell_line, tool) %>%
              unique() %>%
              summarise(total_n = n()) %>% ungroup() %>%
              rename(tool_1 = tool, total_n_1 = total_n)) %>%
  left_join(all_circ %>%
              filter(count_group == 'count ≥ 5') %>%
              group_by(cell_line, tool) %>%
              select(circ_id, cell_line, tool) %>%
              unique() %>%
              summarise(total_n = n()) %>% ungroup() %>%
              rename(tool_2 = tool, total_n_2 = total_n)) %>%
  left_join(all_circ %>%
              filter(count_group == 'count ≥ 5') %>%
              group_by(cell_line, tool) %>%
              select(circ_id, cell_line, tool) %>%
              unique() %>%
              summarise(total_n = n()) %>% ungroup() %>%
              rename(tool_3 = tool, total_n_3 = total_n)) %>%
  # not exactly correct for tools that are combo's of two tools, but will be removed
  mutate(w_val_rate = ((perc_compound_val_1 * total_n_1) + (perc_compound_val_2 * total_n_2) +
                         (perc_compound_val_3 * total_n_3))/(total_n_1 + total_n_2 + total_n_3))

combo_3


#' add tool letters
combo_3 = combo_3 %>%
  left_join(read_tsv('../data/details/tool_annotation.txt') %>%
              select(tool, tool_lt) %>%
              rename(tool_1 = tool, tool_lt_1 = tool_lt)) %>%
  left_join(read_tsv('../data/details/tool_annotation.txt') %>%
              select(tool, tool_lt) %>%
              rename(tool_2 = tool, tool_lt_2 = tool_lt)) %>%
  left_join(read_tsv('../data/details/tool_annotation.txt') %>%
              select(tool, tool_lt) %>%
              rename(tool_3 = tool, tool_lt_3 = tool_lt)) %>%
  mutate(tool_combo = paste(tool_lt_1, tool_lt_2, tool_lt_3, sep = ""),
         tool_combo = map_chr(str_split(tool_combo, ""), ~str_c(str_sort(unique(.x)), collapse = "")))

combo_3

#' remove combo's of two tools
combo_3 = combo_3 %>%
  filter(!nchar(tool_combo) == 2)


#' add total nr of circ for that cell line and calculate perc increase
combo_3 = combo_3 %>%
  left_join(total_cl) %>%
  mutate(perc_increase_t1 = (nr_union - total_n_1)/total_n_1)


#' save as a dataframe for figures
combo_3 %>% write_tsv('../data/Supplementary_Table_8_combo_3tools.txt')
