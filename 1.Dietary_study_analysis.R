# Purpose: To read and process the metabolomics and microbiome data using the maplet pipeline 
# Date: 12-June-2022
# Source of Data : Data received from Dr. David Montrose (Stony Brook university) on 20-March-2022
# Data Description: 17 mice models were fed either a controlled diet (5 mice), 50% reduced protein diet (6 mice), or 50% reduced 
# non-essential amino acid (6 mice) for two weeks. 16s sequencing was performed on fecal samples and metabolomics data performed on feces, liver, and serum.  
# Authors: RKA, NS, KS

rm(list=ls())
library (usethis)
library (devtools)
library(matrixStats)
library (MatrixGenerics)
library(SummarizedExperiment)
library(maplet)
library(tidyverse)
library(magrittr)
library(DT)
library(ggplot2)
library(rstatix)
library(dplyr)
purrr::zap()
file_data <- "Merged_Data.xlsx"
print (file_data )
file_data2 <-  "Merged_Data_2.xlsx"
print (file_data2 )

### Data loading ...
# 1) Loading the data: Fecal 
D_Fecal<- 
  mt_load_xls(file=file_data, sheet="Metabolite_Fecal", samples_in_row=F, id_col="compound") %>%
  # load metabolite (rowData) annotations
  mt_anno_xls(file=file_data, sheet="Metabolomic_Feature",anno_type="features", anno_id_col="compound", data_id_col ="name") %>%
  # load clinical (colData) annotations
  mt_anno_xls(file=file_data, sheet="Samples_Metadata_Fecal", anno_type="samples", anno_id_col ="Sample", data_id_col ="sample") %>%
  # # log assay dimensions and number of columns for both metabolite and clincial annotations
  mt_reporting_data() %>%
  # start timing
  mt_reporting_tic() %>%
  {.}

# 2) Loading the data: Serum 
D_Serum<- 
  mt_load_xls(file=file_data, sheet="Metabolite_Serum", samples_in_row=F, id_col="compound") %>%
  # load metabolite (rowData) annotations
  mt_anno_xls(file=file_data, sheet="Metabolomic_Feature",anno_type="features", anno_id_col="compound", data_id_col ="name") %>%
  # load clinical (colData) annotations
  mt_anno_xls(file=file_data, sheet="Samples_Metadata_Serum", anno_type="samples", anno_id_col ="Sample", data_id_col ="sample") %>%
  # # log assay dimensions and number of columns for both metabolite and clincial annotations
  mt_reporting_data() %>%
  # start timing
  mt_reporting_tic() %>%
  {.}

# 3) Loading the data: Liver
D_Liver<- 
  mt_load_xls(file=file_data, sheet="Metabolite_Liver", samples_in_row=F, id_col="compound") %>%
  # load metabolite (rowData) annotations
  mt_anno_xls(file=file_data, sheet="Metabolomic_Feature",anno_type="features", anno_id_col="compound", data_id_col ="name") %>%
  # load clinical (colData) annotations
  mt_anno_xls(file=file_data, sheet="Sample_Metadata_Liver", anno_type="samples", anno_id_col ="Sample", data_id_col ="sample") %>%
  # # log assay dimensions and number of columns for both metabolite and clinical annotations
  mt_reporting_data() %>%
  # start timing
  mt_reporting_tic() %>%
  {.}

# 4) Loading the data: Microbiome
D_Microbiome <-
  mt_load_xls(file=file_data2, sheet="Species_Count", samples_in_row=F, id_col="species") %>%
  # load (rowData) annotations
  mt_anno_xls(file=file_data2, sheet="Species_Percentage",anno_type="features", anno_id_col="species_P", data_id_col ="name") %>%
  # load clinical (colData) annotations
  mt_anno_xls(file=file_data2, sheet="Microbiome_Metadata", anno_type="samples", anno_id_col ="Sample", data_id_col ="sample") %>%
  # # log assay dimensions and number of columns for both metabolite and clincial annotations
  mt_reporting_data() %>%
  # start timing
  mt_reporting_tic() %>%
  {.}

# Part 1: Analysis of Fecal Data set 

# Part 1.1: MISSINGNESS ANALYSIS 

D_Fecal <- D_Fecal %>%
  mt_pre_zero_to_na() %>%
  mt_reporting_heading(heading = "Missingness analysis", lvl = 1) %>%
  # section text
  mt_reporting_text(text = "Perform missingness analysis to determine if NAs significantly accumulate in one of the treatmen
                    groups. Adjust output of test using multiple testing correction.") %>%
  # compute Fisher's exact test
  mt_stats_univ_missingness(in_col="Group", stat_name="missingness") %>%
  # create p-value qq plot
  mt_plots_pval_qq(stat_name = "missingness") %>%
  # apply multiple testing correction
  #   alternative function: mt_post_multtest_effdim
  mt_post_multtest(stat_name="missingness", method="BH") %>%
  {.}


# Part 1.2: FILTERING MISSING VALUES

D_Fecal <- D_Fecal %>%
  mt_modify_filter_samples(filter = !is.na(Group)) %>%
  mt_anno_apply(anno_type = "samples", col_name = "Group", fun = as.factor) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  mt_pre_filter_missingness(feat_max = 0.2) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  {.}


# Part 1.3: PREPROCESSING: NORMALIZATION

D_Fecal <- D_Fecal %>%
  # heading for html file
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot sample boxplots before normalization, perform quotient
                    normalization, plot boxplot with dilution factors from quotient normalization, plot sample boxplot after
                    normalization, log transform the data, impute missing data using min value, plot sample boxplot after imputation,
                    detect outliers, log dataset info, write pre-processed data to file.") %>%
  mt_plots_sample_boxplot(color=Group, title = "Original", plot_logged = T) %>%
  mt_pre_norm_quot(feat_max = 0.2, ref_samples = Group=="Control") %>%
  mt_plots_dilution_factor(in_col="Group") %>%
  # plot sample boxplots after normalization
  mt_plots_sample_boxplot(color=Group, title = "After normalization", plot_logged = T) %>%
  # log transform
  mt_pre_trans_log() %>%
  # impute missing values using min
  mt_pre_impute_min() %>%
  # plot sample boxplot after imputation
  mt_plots_sample_boxplot(color=Group, title = "After imputation", plot_logged = T) %>%
  mt_pre_outlier_detection_univariate() %>%
  mt_reporting_data() %>%
  
  {.}

# Part 1.4: GLOBAL STATISTICS 

assay(D_Fecal)[ which(is.na(assay(D_Fecal)), arr.ind = TRUE) ] = 0

D_Fecal <- D_Fecal %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  
  # plot PCA tissue Group
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", color= Group, exp_var_plot=T, size=2.5, ggadd= theme(panel.background = element_rect(fill = "white", colour = "grey50"))) %>%
  
  {.}



D_Fecal <- D_Fecal %>%
  # heading for html file
  mt_reporting_heading(heading = "D_FecalGroup analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Group,
                   stat_name = "D_FecalGroup met") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "D_FecalGroup met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "D_FecalGroup met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "D_FecalGroup met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "D_FecalGroupAnalysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "D_FecalGroup met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "D_FecalGroup met",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="D_FecalGroup met",
                       x                  = Group,
                       fill               = Group,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
  {.}

# Part 2: Analysis of Serum Data set 

# Part 2.1: MISSINGNESS ANALYSIS 

D_Serum <- D_Serum %>%
  mt_pre_zero_to_na() %>%
  mt_reporting_heading(heading = "Missingness analysis", lvl = 1) %>%
  # section text
  mt_reporting_text(text = "Perform missingness analysis to determine if NAs significantly accumulate in one of the treatmen
                    groups. Adjust output of test using multiple testing correction.") %>%
  # compute Fisher's exact test
  mt_stats_univ_missingness(in_col="Group", stat_name="missingness") %>%
  # create p-value qq plot
  mt_plots_pval_qq(stat_name = "missingness") %>%
  # apply multiple testing correction
  #   alternative function: mt_post_multtest_effdim
  mt_post_multtest(stat_name="missingness", method="BH") %>%
  {.}


# Part 2.2: FILTERING MISSING VALUES 

D_Serum <- D_Serum %>%
  mt_modify_filter_samples(filter = !is.na(Group)) %>%
  mt_anno_apply(anno_type = "samples", col_name = "Group", fun = as.factor) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  mt_pre_filter_missingness(feat_max = 0.2) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  {.}


# Part 2.3: PREPROCESSING: NORMALIZATION 

D_Serum <- D_Serum %>%
  # heading for html file
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot sample boxplots before normalization, perform quotient
                    normalization, plot boxplot with dilution factors from quotient normalization, plot sample boxplot after
                    normalization, log transform the data, impute missing data using min value, plot sample boxplot after imputation,
                    detect outliers, log dataset info, write pre-processed data to file.") %>%
  mt_plots_sample_boxplot(color=Group, title = "Original", plot_logged = T) %>%
  mt_pre_norm_quot(feat_max = 0.2, ref_samples = Group=="Control") %>%
  mt_plots_dilution_factor(in_col="Group") %>%
  # plot sample boxplots after normalization
  mt_plots_sample_boxplot(color=Group, title = "After normalization", plot_logged = T) %>%
  # log transform
  mt_pre_trans_log() %>%
  # impute missing values using min
  mt_pre_impute_min() %>%
  # plot sample boxplot after imputation
  mt_plots_sample_boxplot(color=Group, title = "After imputation", plot_logged = T) %>%
  mt_pre_outlier_detection_univariate() %>%
  mt_reporting_data() %>%
  
  {.}


# Part 2.4: GLOBAL STATISTICS 

assay(D_Serum)[ which(is.na(assay(D_Serum)), arr.ind = TRUE) ] = 0

D_Serum <- D_Serum %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  
  # plot PCA tissue Group
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", color= Group, exp_var_plot=T, size=2.5, ggadd= theme(panel.background = element_rect(fill = "white", colour = "grey50"))) %>%
  {.}


D_Serum <- D_Serum %>% 
  # heading for html file
  mt_reporting_heading(heading = "D_SerumGroup analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Group,
                   stat_name = "D_SerumGroup met") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "D_SerumGroup met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "D_SerumGroup met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "D_SerumGroup met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "D_SerumGroupAnalysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "D_SerumGroup met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "D_SerumGroup met",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="D_SerumGroup met",
                       x                  = Group,
                       fill               = Group,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
  {.}


# Part 3: Analysis of Liver Data set 

# Part 3.1: MISSINGNESS ANALYSIS 

D_Liver <- D_Liver %>%
  mt_pre_zero_to_na() %>%
  mt_reporting_heading(heading = "Missingness analysis", lvl = 1) %>%
  # section text
  mt_reporting_text(text = "Perform missingness analysis to determine if NAs significantly accumulate in one of the treatmen
                    groups. Adjust output of test using multiple testing correction.") %>%
  # compute Fisher's exact test
  mt_stats_univ_missingness(in_col="Group", stat_name="missingness") %>%
  # create p-value qq plot
  mt_plots_pval_qq(stat_name = "missingness") %>%
  # apply multiple testing correction
  #   alternative function: mt_post_multtest_effdim
  mt_post_multtest(stat_name="missingness", method="BH") %>%
  {.}


# Part 3.2: FILTERING MISSING VALUES  

D_Liver <- D_Liver %>%
  mt_modify_filter_samples(filter = !is.na(Group)) %>%
  mt_anno_apply(anno_type = "samples", col_name = "Group", fun = as.factor) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  mt_pre_filter_missingness(feat_max = 0.2) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  {.}


# Part 3.3: PREPROCESSING: NORMALIZATION 

D_Liver <- D_Liver %>%
  # heading for html file
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot sample boxplots before normalization, perform quotient
                    normalization, plot boxplot with dilution factors from quotient normalization, plot sample boxplot after
                    normalization, log transform the data, impute missing data using min value, plot sample boxplot after imputation,
                    detect outliers, log dataset info, write pre-processed data to file.") %>%
  mt_plots_sample_boxplot(color=Group, title = "Original", plot_logged = T) %>%
  mt_pre_norm_quot(feat_max = 0.2, ref_samples = Group=="Control") %>%
  mt_plots_dilution_factor(in_col="Group") %>%
  # plot sample boxplots after normalization
  mt_plots_sample_boxplot(color=Group, title = "After normalization", plot_logged = T) %>%
  # log transform
  mt_pre_trans_log() %>%
  # impute missing values using min
  mt_pre_impute_min() %>%
  # plot sample boxplot after imputation
  mt_plots_sample_boxplot(color=Group, title = "After imputation", plot_logged = T) %>%
  mt_pre_outlier_detection_univariate() %>%
  mt_reporting_data() %>%
  
  {.}


# Part 3.4: GLOBAL STATISTICS 

assay(D_Liver)[ which(is.na(assay(D_Liver)), arr.ind = TRUE) ] = 0

D_Liver <- D_Liver %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  
  # plot PCA tissue Group
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", color= Group, exp_var_plot=T, size=2.5, ggadd= theme(panel.background = element_rect(fill = "white", colour = "grey50"))) %>%
  {.}


D_Liver <- D_Liver %>%
  # heading for html file
  mt_reporting_heading(heading = "D_LiverGroup analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Group,
                   stat_name = "D_LiverGroup met") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "D_LiverGroup met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "D_LiverGroup met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "D_LiverGroup met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "D_LiverGroupAnalysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "D_LiverGroup met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "D_LiverGroup met",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="D_LiverGroup met",
                       x                  = Group,
                       fill               = Group,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
  {.}

# Part 4: Analysis of Microbiome Data set 

# Part 4.1: MISSINGNESS ANALYSIS 

D_Microbiome  <- D_Microbiome  %>%
  mt_pre_zero_to_na() %>%
  mt_reporting_heading(heading = "Missingness analysis", lvl = 1) %>%
  # section text
  mt_reporting_text(text = "Perform missingness analysis to determine if NAs significantly accumulate in one of the treatmen
                    groups. Adjust output of test using multiple testing correction.") %>%
  # compute Fisher's exact test
  mt_stats_univ_missingness(in_col="Group", stat_name="missingness") %>%
  # create p-value qq plot
  mt_plots_pval_qq(stat_name = "missingness") %>%
  # apply multiple testing correction
  #   alternative function: mt_post_multtest_effdim
  mt_post_multtest(stat_name="missingness", method="BH") %>%
  {.}


# Part 4.2: FILTERING MISSING VALUES

D_Microbiome  <- D_Microbiome  %>%
  mt_modify_filter_samples(filter = !is.na(Group)) %>%
  mt_anno_apply(anno_type = "samples", col_name = "Group", fun = as.factor) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  mt_pre_filter_missingness(feat_max = 0.2) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  {.}

# from 95 to 64 

# Part 4.3: PREPROCESSING: NORMALIZATION

D_Microbiome  <- D_Microbiome  %>%
  # heading for html file
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot sample boxplots before normalization, perform quotient
                    normalization, plot boxplot with dilution factors from quotient normalization, plot sample boxplot after
                    normalization, log transform the data, impute missing data using min value, plot sample boxplot after imputation,
                    detect outliers, log dataset info, write pre-processed data to file.") %>%
  mt_plots_sample_boxplot(color=Group, title = "Original", plot_logged = T) %>%
  mt_pre_norm_quot(feat_max = 0.2, ref_samples = Group=="Control") %>%
  mt_plots_dilution_factor(in_col="Group") %>%
  # plot sample boxplots after normalization
  mt_plots_sample_boxplot(color=Group, title = "After normalization", plot_logged = T) %>%
  # log transform
  mt_pre_trans_log() %>%
  # impute missing values using min
  mt_pre_impute_min() %>%
  # plot sample boxplot after imputation
  mt_plots_sample_boxplot(color=Group, title = "After imputation", plot_logged = T) %>%
  mt_pre_outlier_detection_univariate() %>%
  mt_reporting_data() %>%
  
  {.}

# Part 4.4: GLOBAL STATISTICS 

assay(D_Microbiome)[ which(is.na(assay(D_Microbiome)), arr.ind = TRUE) ] = 0

D_Microbiome  <- D_Microbiome  %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  
  # plot PCA tissue Group
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", color= Group, exp_var_plot=T, size=2.5, ggadd= theme(panel.background = element_rect(fill = "white", colour = "grey50"))) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", color=Group, exp_var_plot=T, size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", pc1=1, pc2=3, color=Group,size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", pc1=2, pc2=3, color=Group,size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", pc1=4, pc2=5, color=Group,size=2.5, ggadd=scale_size_identity()) %>%
  {.}


D_Microbiome  <- D_Microbiome  %>%
  # heading for html file
  mt_reporting_heading(heading = "D_MGroup analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Group,
                   stat_name = "D_MGroup analysis") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "D_MGroup  met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "D_MGroup analysis", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "D_MGroup analysis", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "D_MGroup analysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "D_MGroup analysis") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "D_MGroup analysis",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="D_MGroup analysis",
                       x                  = Group,
                       fill               = Group,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
  {.}


# Part 5: Merge all datasets together 

# merging functions: 

rmerge_SE <- function(D1, D2) {
  
  cat("merging two SE by row:\n")
  cat("D1:\n")
  print(D1)
  cat("D2:\n")
  print(D2)
  
  # check if the colnames match
  # a future version could use left_join and pass through its options
  if (! sum(colnames(D1) == colnames(D2)) == dim(D1)[2]) {
    stop("colnames(D1) != colnames(D2)")
  }
  
  # remove all colData that does not have the same variable name and class and same content for each sample
  # if data contains NA's, it will be removed ... this should be improved in a future version
  # but I knew of no quick way of comparing vectors for identity if they contain NA's
  # (there should be a simple routine for this somewhere),
  # sum(cdata1 == cdata2, na.rm = TRUE) != dim(D1)[2] is not very elegant
  # a future version could have the option to keep non-matching variables and rename them 
  # to .1 and .2 and fill the other data with NA's
  # if the two SE's were read by Maplet, then also the merging of the meta information needs to be treated
  # i.e. for now it reports only the first loaded file in the HTML output
  for (i in intersect(names(colData(D1)),names(colData(D2))) ) {
    # cat("working on", i, "\n")
    cdata1 = D1[[i]]
    cdata2 = D2[[i]]
    lidentical = TRUE
    if (class(cdata1) != class(cdata2)) {
      lidentical = FALSE
      # cat("not identical - different class type\n")
    } else if (sum(cdata1 == cdata2, na.rm = TRUE) != dim(D1)[2]) {
      lidentical = FALSE
      # cat("not identical - different content\n")
    }
    if (! lidentical) {
      cat("deleting variable with non-matching values in D1 and D2:", i, "\n")
      colData(D1)[which(names(colData(D1)) == i)] = NULL
      colData(D2)[which(names(colData(D2)) == i)] = NULL
    }
    
  }
  
  # remove variables that are only in one data set, 
  # a future version could keep them and fill in NA's
  for (i in setdiff(names(colData(D1)),names(colData(D2))) ) {
    cat("deleting variable that is only in D1:", i, "\n")
    colData(D1)[which(names(colData(D1)) == i)] = NULL
  }  
  for (i in setdiff(names(colData(D2)),names(colData(D1))) ) {
    cat("deleting variable that is only in D2:", i, "\n")
    colData(D2)[which(names(colData(D2)) == i)] = NULL
  }  
  
  
  # treat row data: remove variables that are in only one data set
  # future: could be filled with NA's
  # also make sure they are of the same class
  
  for (i in setdiff(names(rowData(D1)),names(rowData(D2))) ) {
    cat("deleting row variable that is only in D1:", i, "\n")
    rowData(D1)[which(names(rowData(D1)) == i)] = NULL
  }
  
  for (i in setdiff(names(rowData(D2)),names(rowData(D1))) ) {
    cat("deleting row variable that is only in D2:", i, "\n")
    rowData(D2)[which(names(rowData(D2)) == i)] = NULL
  }  
  
  for (i in intersect(names(rowData(D1)),names(rowData(D2))) ) {
    # cat("working on", i, "\n")
    if (class(D1[[i]]) != class(D1[[i]])) {
      cat("deleting row variable with different class in D1 and D2:", i, "\n")
      rowData(D1)[which(names(rowData(D1)) == i)] = NULL
      rowData(D2)[which(names(rowData(D2)) == i)] = NULL
    }  
  }  
  
  # add a variable to id datasets
  rowData(D1)$original_dataset = rep("D1", dim(D1)[1])
  rowData(D2)$original_dataset = rep("D2", dim(D2)[1])
  
  # return merged SE
  rbind(D1,D2)
  
}

cmerge_SE <- function(D1, D2) {
  
  # test code
  # D1 = D_muscle
  # D2 = D_fat
  # cmerge_SE(D_muscle, D_fat)
  
  cat("merging two SE by column:\n")
  cat("D1:\n")
  print(D1)
  cat("D2:\n")
  print(D2)
  
  # check if the colnames match
  # a future version could use left_join and pass through its options
  if (! sum(rownames(D1) == rownames(D2)) == dim(D1)[1]) {
    stop("rownames(D1) != rownames(D2)")
  }
  
  # remove all rowData that does not have the same variable name and class and same content for each sample
  for (i in intersect(names(rowData(D1)),names(rowData(D2))) ) {
    # cat("working on", i, "\n")
    cdata1 = rowData(D1)[[i]]
    cdata2 = rowData(D2)[[i]]
    lidentical = TRUE
    if (class(cdata1) != class(cdata2)) {
      lidentical = FALSE
      # cat("not identical - different class type\n")
    } else if (sum(cdata1 == cdata2, na.rm = TRUE) != dim(D1)[1]) {
      lidentical = FALSE
      # cat("not identical - different content\n")
    }
    if (! lidentical) {
      cat("deleting variable with non-matching values in D1 and D2:", i, "\n")
      rowData(D1)[which(names(rowData(D1)) == i)] = NULL
      rowData(D2)[which(names(rowData(D2)) == i)] = NULL
    }
    
  }
  
  # remove variables that are only in one data set, 
  # a future version could keep them and fill in NA's
  for (i in setdiff(names(rowData(D1)),names(rowData(D2))) ) {
    cat("deleting variable that is only in D1:", i, "\n")
    rowData(D1)[which(names(rowData(D1)) == i)] = NULL
  }  
  for (i in setdiff(names(rowData(D2)),names(rowData(D1))) ) {
    cat("deleting variable that is only in D2:", i, "\n")
    rowData(D2)[which(names(rowData(D2)) == i)] = NULL
  }  
  
  # treat row data: remove variables that are in only one data set
  # future: could be filled with NA's
  # also make sure they are of the same class
  
  for (i in setdiff(names(colData(D1)),names(colData(D2))) ) {
    cat("deleting col variable that is only in D1:", i, "\n")
    colData(D1)[which(names(colData(D1)) == i)] = NULL
  }
  
  for (i in setdiff(names(colData(D2)),names(colData(D1))) ) {
    cat("deleting col variable that is only in D2:", i, "\n")
    colData(D2)[which(names(colData(D2)) == i)] = NULL
  }  
  
  for (i in intersect(names(colData(D1)),names(colData(D2))) ) {
    # cat("working on", i, "\n")
    if (class(colData(D1)[[i]]) != class(colData(D2)[[i]])) {
      cat("deleting col variable with different class in D1 and D2:", i, "\n")
      colData(D1)[which(names(colData(D1)) == i)] = NULL
      colData(D2)[which(names(colData(D2)) == i)] = NULL
    }  
  }  
  
  # add a variable to id datasets
  colData(D1)$original_dataset = rep("D1", dim(D1)[2])
  colData(D2)$original_dataset = rep("D2", dim(D2)[2])
  
  # return merged SE
  cbind(D1,D2)
  
}

# ColData merge for Liver, Serum, and Fecal Matrices: 

# 5.1 Pre- merging steps:

rowData(D_Fecal)$HMDB[which(is.na(rowData(D_Fecal)$HMDB))]= "Non"
rowData(D_Serum)$HMDB[which(is.na(rowData(D_Serum)$HMDB))]= "Non"
rowData(D_Liver)$HMDB[which(is.na(rowData(D_Liver)$HMDB))]= "Non"
D_Fecal2<- D_Fecal%>% mt_pre_trans_scale()
D_Serum2<- D_Serum%>% mt_pre_trans_scale()
D_Liver2<- D_Liver%>% mt_pre_trans_scale()
colnames(D_Fecal2)=paste(colnames(D_Fecal2),"Fecal")
colnames(D_Serum2)=paste(colnames(D_Serum2),"Serum")
colnames(D_Liver2)=paste(colnames(D_Liver2),"Liver")
shared= intersect(intersect(rownames(D_Fecal2),rownames(D_Liver2)),rownames(D_Serum2))
D_Fecal2= D_Fecal2[left_join(data.frame(id=shared),data.frame(id=rownames(D_Fecal2),ix=seq(dim(D_Fecal2)[1])))$ix,]
D_Liver2= D_Liver2[left_join(data.frame(id=shared),data.frame(id=rownames(D_Liver2),ix=seq(dim(D_Liver2)[1])))$ix,]
D_Serum2= D_Serum2[left_join(data.frame(id=shared),data.frame(id=rownames(D_Serum2),ix=seq(dim(D_Serum2)[1])))$ix,]

# 5.2 Merge the three matrices: 

D_merged = cmerge_SE(D_Fecal2, D_Serum2)
D_merged2 = cmerge_SE(D_merged, D_Liver2)
D_merged2$tissue=sub(".* ", "",colnames(D_merged2))
colnames(D_merged2)=make.names(colnames(D_merged2))
metadata(D_merged2)= list()

D <- D_merged2 %>%
  mt_clean_validate_se()%>%
  mt_clean_remove_results()%>%
  {.}


# RowData merge: 

D_Fecal2<- D_Fecal
D_Serum2<- D_Serum
D_Liver2<- D_Liver
#intersect (colnames(D_Fecal2),colnames(D_Liver2))
D_Fecal2<- D_Fecal2[,which(colnames(D_Fecal2)!="CTR_2")]
D_Serum2<- D_Serum2[,which(colnames(D_Serum2)!="CTR_2")]
rownames(D_Fecal2)=paste(rownames(D_Fecal2),"Fecal")
rownames(D_Serum2)=paste(rownames(D_Serum2),"Serum")
rownames(D_Liver2)=paste(rownames(D_Liver2),"Liver")
D_rmerged = rmerge_SE(D_Fecal2, D_Serum2)
D_rmerged2 = rmerge_SE(D_rmerged, D_Liver2)
rowData(D_rmerged2)$tissue=sub(".* ", "",rownames(D_rmerged2))
rownames(D_rmerged2)=make.names(rownames(D_rmerged2))
head(rowData(D_rmerged2)$HMDB)
names((rowData(D_rmerged2)))
head(rowData(D_rmerged2)$name)
head(rowData(D_rmerged2)$compound)
rowData(D_rmerged2)$name= paste(rowData(D_rmerged2)$name,rowData(D_rmerged2)$tissue)
colnames(D_Fecal2)==colnames(D_Serum2)

# Microbiome RowData merge: 

D_Microbiome<- D_Microbiome[,which(colnames(D_Microbiome)!="CTR_2")]
rownames(D_Microbiome)=paste(rownames(D_Microbiome),"Microbiome")
D_rmerged3 = rmerge_SE(D_rmerged2, D_Microbiome)
rownames(D_rmerged3)=make.names(rownames(D_rmerged3))
rownames(D_rmerged3)


# Part 6: Analysis of all Data set 

# ColData: 
#assay(D)[ which(is.na(assay(D)), arr.ind = TRUE) ] = 0

D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  # plot PCA tissue Group
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", color=Group, shape = tissue, exp_var_plot=T, size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", pc1=1, pc2=3, color=Group, shape = tissue,size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", pc1=2, pc2=3, color=Group, shape = tissue,size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Group", pc1=4, pc2=5, color=Group, shape = tissue,size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_heatmap(scale_data = T, annotation_col = c("Group"), 
                   clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
  {.}


 

D1<- D
D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "D_merged2 analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~ Group,
                   stat_name = "D_merged2 met") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "D_merged2 met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "D_merged2 met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "D_merged2 met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "D_merged2Analysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "D_merged2 met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "D_merged2 met",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="D_merged2 met",
                       x                  = Group,
                       fill               = Group,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
   {.}


# RowData : 

metadata(D_rmerged2)= list()
D_rmerged2 <- D_rmerged2 %>%
  mt_clean_validate_se()%>%
  mt_clean_remove_results()%>%
  {.}

D_rmerged2 <- D_rmerged2 %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  # plot PCA 
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Rownames merged", color= Group, exp_var_plot=T, size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Rownames merged", pc1=1, pc2=3, color=Group, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Rownames merged", pc1=2, pc2=3, color=Group, ggadd=scale_size_identity()) %>%
  mt_plots_heatmap(scale_data = T, annotation_col = c("Group"), 
                   clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
  {.}

# correlation analysis: 

assay(D_rmerged2)
x=cor(t(assay(D_rmerged2)))
str(x)
hist(x)
#r2=x^2
#hist(x^2)
y= reshape2::melt(x)
y$Var1=as.character(y$Var1)
y$Var2=as.character(y$Var2)
ix= which(y$Var1<y$Var2)
y=y[ix,]
y$r2=y$value^2
ix2=sort(y$r2,decreasing= TRUE,index.return=TRUE)$ix
y=y[ix2,]
head(y,1000)
write.table(y,file="correlation.tsv",sep="\t",row.names=FALSE)

# correlation_microbiome_only
assay(D_Microbiome)
x=cor(t(assay(D_Microbiome)))
str(x)
hist(x)
#r2=x^2
#hist(x^2)
y= reshape2::melt(x)
y$Var1=as.character(y$Var1)
y$Var2=as.character(y$Var2)
ix= which(y$Var1<y$Var2)
y=y[ix,]
y$r2=y$value^2
ix2=sort(y$r2,decreasing= TRUE,index.return=TRUE)$ix
y=y[ix2,]
head(y,1000)
write.table(y,file="correlationMOnly.tsv",sep="\t",row.names=FALSE)

# GGM model: 
# D_rmerged2 <- D_rmerged2 %>%
#   mt_stats_cormat_genenet(stat_name="pcor_row") %>%
#   mt_plots_net(stat_name = "pcor_row") %>%
#   {.}



D_rmerged2 <- D_rmerged2 %>%
  # heading for html file
  mt_reporting_heading(heading = "D_rmerged2 analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~Group ,
                   stat_name = "D_rmerged2 met") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "D_rmerged2 met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "D_rmerged2 met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "D_rmerged2 met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "D_rmerged2Analysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "D_rmerged2 met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "D_rmerged2 met",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  # mt_plots_box_scatter(stat_name          ="D_rmerged2 met",
  #                      x                  = Group,
  #                      fill               = Group,
  #                      plot_type          = "box",
  #                      feat_filter       = p.adj < 0.05,
  #                      feat_sort         = p.value,
  #                      annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
  {.}

# creating additional variables to compare between the administered diet

table(D_rmerged3$Group)
str(D_rmerged3$Group)

# group control vs reduced protein 
D_rmerged3$Gctr_prot=rep(NA,dim(D_rmerged3)[2])
D_rmerged3$Gctr_prot[which(D_rmerged3$Group=="Control")]="control"
D_rmerged3$Gctr_prot[which(D_rmerged3$Group=="Reduced Protein")]="Reduced Protein"
table(D_rmerged3$Gctr_prot,useNA="ifany")

# group control vs NEAA
D_rmerged3$Gctr_NEAA=rep(NA,dim(D_rmerged3)[2])
D_rmerged3$Gctr_NEAA[which(D_rmerged3$Group=="Control")]="control"
D_rmerged3$Gctr_NEAA[which(D_rmerged3$Group=="Reduced NEAA")]="Reduced NEAA"
table(D_rmerged3$Gctr_NEAA,useNA="ifany")

# group control vs both diet
D_rmerged3$Gctr_both=rep(NA,dim(D_rmerged3)[2])
D_rmerged3$Gctr_both[which(D_rmerged3$Group=="Control")]="control"
D_rmerged3$Gctr_both[which(D_rmerged3$Group!="Control")]="both"
table(D_rmerged3$Gctr_both,useNA="ifany")

# reduced protein vs. NEAA 
D_rmerged3$prot_NEAA=rep(NA,dim(D_rmerged3)[2])
D_rmerged3$prot_NEAA[which(D_rmerged3$Group=="Reduced Protein")]="Reduced Protein"
D_rmerged3$prot_NEAA[which(D_rmerged3$Group=="Reduced NEAA")]="Reduced NEAA"

# STATISTICAL ANALYSIS: 

D_rmerged_ALL <- D_rmerged3 %>%
  mt_clean_remove_results()%>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~Group ,
                   stat_name = "Group") %>%
  mt_stats_univ_lm(formula = ~Gctr_NEAA ,
                   stat_name = "Gctr_NEAA") %>%
  mt_stats_univ_lm(formula = ~Gctr_prot ,
                   stat_name = "Gctr_prot") %>%
  mt_stats_univ_lm(formula = ~Gctr_both ,
                   stat_name = "Gctr_both") %>%
  mt_stats_univ_lm(formula = ~prot_NEAA ,
                   stat_name = "prot_NEAA") %>%
  mt_write_stats(file = "D_rmerged_ALL_Analysis.xlsx",stat_list=c("Group","Gctr_NEAA","Gctr_prot","Gctr_both","prot_NEAA")) %>%
  {.}

# GGM model: 
# D_rmerged_ALL <- D_rmerged_ALL %>%
#   mt_stats_cormat_genenet(stat_name="pcor_row_ALL") %>%
#   mt_plots_net(stat_name = "pcor_row_ALL") %>%
#   {.}

# Microbiome merged analysis 

metadata(D_rmerged3)= list()
D_rmerged3 <- D_rmerged3 %>%
  mt_clean_validate_se()%>%
  mt_clean_remove_results()%>%
  {.}

D_rmerged3 <- D_rmerged3 %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  # plot PCA 
  mt_plots_pca(scale_data = T, title = "Scaled PCA - Microbiome Rownames merged", color= Group, exp_var_plot=T, size=0.5, ggadd= theme(panel.background = element_rect(fill = "white", colour = "grey50"))) %>%
  mt_plots_heatmap(scale_data = T, annotation_col = c("Group"), 
                   clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
  {.}
rowData(D_rmerged3)$name

# correlation analysis: 

assay(D_rmerged3)
x=cor(t(assay(D_rmerged3)))
str(x)
hist(x)
#r2=x^2
#hist(x^2)
y= reshape2::melt(x)
y$Var1=as.character(y$Var1)
y$Var2=as.character(y$Var2)
ix= which(y$Var1<y$Var2)
y=y[ix,]
y$r2=y$value^2
ix2=sort(y$r2,decreasing= TRUE,index.return=TRUE)$ix
y=y[ix2,]
head(y,1000)
write.table(y,file="correlation_microbiome.tsv",sep="\t",row.names=FALSE)


# GGM model: 
# D_rmerged3 <- D_rmerged3 %>%
#   mt_stats_cormat_genenet(stat_name="pcor_row") %>%
#   mt_plots_net(stat_name = "pcor_row") %>%
#   {.}


# STATISTICAL ANALYSIS, OUTCOME: Group, METHOD: LINEAR REGRESSION (t-test) 
D_rmerged3 <- D_rmerged3 %>%
  # heading for html file
  mt_reporting_heading(heading = "D_rmerged3 analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  mt_stats_univ_lm(formula = ~Group ,
                   stat_name = "D_rmerged3 met") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "D_rmerged3 met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "D_rmerged3 met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "D_rmerged3 met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "D_rmerged3Analysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "D_rmerged3 met") %>%
  # volcano plot as overview of results
  mt_plots_volcano(stat_name = "D_rmerged3 met",
                   x = statistic,
                   feat_filter = p.adj < 0.05,
                   colour       = p.adj < 0.05) %>%
  # boxplot
  # mt_plots_box_scatter(stat_name          ="D_rmerged3 met",
  #                      x                  = Group,
  #                      fill               = Group,
  #                      plot_type          = "box",
  #                      feat_filter       = p.adj < 0.05,
  #                      feat_sort         = p.value,
  #                      annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
   {.}



# Creating the reports for the tested SE

D_Serum <- D_Serum %>% mt_reporting_html(file = "D_SerumAnalysis.html",
                                         title = "MiceData - Statistical Analysis")

D_Fecal <- D_Fecal %>% mt_reporting_html(file = "D_FecalAnalysis.html",
                                         title = "MiceData - Statistical Analysis")

D_Liver <- D_Liver %>% mt_reporting_html(file = "D_LiverAnalysis.html",
                                         title = "MiceData - Statistical Analysis")

D_Microbiome <- D_Microbiome %>% mt_reporting_html(file = "D_MicrobiomeAnalysis.html",
                                         title = "MiceData - Statistical Analysis")

D1<- D1 %>% mt_reporting_html(file = "PCA_ColData_Merged_SE.html",
                              title = "Example Pipeline Metabolome/Microbiome - Statistical Analysis")

D_rmerged2 <- D_rmerged2 %>% mt_reporting_html(file = "PCA_Rowdata_Merged_SE.html",
                                               title = "Example Pipeline Metabolome/Microbiome - Statistical Analysis")

D_rmerged3 <- D_rmerged3 %>% mt_reporting_html(file = "PCA_Microbiome_Rowdata_Merged_SE.html",
                                               title = "Example Pipeline Metabolome/Microbiome - Statistical Analysis")

D_rmerged4 <- D_rmerged4 %>% mt_reporting_html(file = "Rowdata_Merged_SE.html",
                                               title = "Example Pipeline Metabolome/Microbiome - Statistical Analysis")

