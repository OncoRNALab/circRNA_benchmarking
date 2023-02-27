# circRNA_benchmarking

This repository contains all data and scripts used to generate the numbers and figures for the circRNA detection tool benchmarking paper (currently on [BioRxiv](https://www.biorxiv.org/content/10.1101/2022.12.06.519083v1)).

The `data` folder contains
- the `Supplementary_Table_2_all_circRNAs.txt` file, which contains all predicted circRNAs (1.1 million) in this study with their annotation.
- the `Supplementary_Table_3_selected_circRNAs.txt` file, which contains all 1560 circRNAs selected for validation, with their initial detection information (tool, BSJ count), their primer information (including FWD and REV primer sequence), results from three validation methods (Cq value with and without RNase R, Cq difference, and amplicon sequencing percent on-target amplification), validation metrics, and annotation information. This file was generated by the `01_calculate_val_rates.R` script.
- the `Supplementary_Table_4_precision_values.txt` files, which contains the validation metrics (per-methods precision, compound precision, theoretical number of TP circRNAs, estimated sensitivity) for each tool. This file was generated using the `01_calculate_val_rates.R` script.
- the  `Supplementary_Table_5_combo_2tools.txt` file, which contains the number of circRNAs in the intersection and union of each combination of two tools, per cell line sample (only for circRNAs with BSJ count ≥ 5). This file was generated by the `04_combination_tools.R` script.
- the `Supplementary_Table_6_combo_3tools.txt` file, which contains the number of circRNAs in the intersection and union of each combination of three tools, per cell line sample (only for circRNAs with BSJ count ≥ 5). This file was generated by the `04_combination_tools.R` script.
- the `details` folder, which contains some files needed for the following scripts to generate some of the Supplementary Figures and Tables. `circ_db_hg38.txt` is a table with all circRNAs in all circRNA databases from a previous [publication](https://academic.oup.com/bib/article/22/1/288/5717788).


The `data_analysis` folder contains
- `01_calculate_val_rates.R` file, which contains the calculations of the validation metrics (per-methods precision, compound precision, theoretical number of TP circRNAs, estimated sensitivity) and generates `Supplementary_Table_3_selected_circRNAs.txt` and `Supplementary_Table_4_precision_values.txt`.
- `02_calculations_paper.R` file, which contains all calculations reported in the manuscript.
- `03_annotation_and_validation.R` file, which contains all calulcations described in the paragraph *Comparing precision values in function of circRNA annotation.*
- `04_combination_tools.R` file, which contains all calculation for the union and intersection of two or three tools and generates `Supplementary_Table_5_combo_2tools.txt` and `Supplementary_Table_6_combo_3tools.txt`.

The `figure_generating` folder contains the R scripts and R markdowns to generate all Figures and Supplementary Figures in the manuscript.


## citation
Large-scale benchmarking of circRNA detection tools reveals large differences in sensitivity but not in precision
Marieke Vromman, Jasper Anckaert, Stefania Bortoluzzi, Alessia Buratin, Chia-Ying Chen, Qinjie Chu, Trees-Juen Chuang, Roozbeh Dehghannasiri, Christoph Dieterich, Xin Dong, Paul Flicek, Enrico Gaffo, Wanjun Gu, Chunjiang He, Steve Hoffmann, Osagie Izuogu, Michael S. Jackson, Tobias Jakobi, Eric C. Lai, Justine Nuytens, Julia Salzman, Mauro Santibanez-Koref, Peter Stadler, Olivier Thas, Eveline Vanden Eynde, Kimberly Verniers, Guoxia Wen, Jakub Westholm, Li Yang, Chu-Yu Ye, Nurten Yigit, Guo-Hua Yuan, Jinyang Zhang, Fangqing Zhao, Jo Vandesompele, Pieter-Jan Volders
bioRxiv 2022.12.06.519083; doi: https://doi.org/10.1101/2022.12.06.519083
