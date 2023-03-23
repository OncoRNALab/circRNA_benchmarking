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

#' #' clean up dataframes and add validation rates
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

#' #' add total nr of circ per tool
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

#' #' add tool letters
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


#' save as a dataframe for Supplementary table 5

combo %>% write_tsv('../data/Supplementary_Table_7_combo_2tools.txt')


#' ## Generate a top 10 list of combinations (Supplementary Table 9)

#' select only one combination when there are doubles + remove the combinations of twice the same tool
combo_top = combo %>%
  filter(!tool_1 == tool_2) %>%
  group_by(tool_combo, cell_line) %>%
  slice(1) %>%
  ungroup() 

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
  unique() %>%
  select("tool_1", "tool_2", "nr_union_HLF", "nr_union_HLF_rank", "nr_union_NCI_H23", 
         "nr_union_NCI_H23_rank", "nr_union_SW480", "nr_union_SW480_rank" , "w_precision_HLF", 
         "w_precision_HLF_rank",  "w_precision_NCI_H23", "w_precision_NCI_H23_rank",
         "w_precision_SW480", "w_precision_SW480_rank")

#' save as sup table 9
combo_top %>% write_tsv('../data/Supplementary_Table_9_top_tool_combinations.txt')


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

#' save as a dataframe for figures
combo_3 %>% write_tsv('../data/Supplementary_Table_8_combo_3tools.txt')
