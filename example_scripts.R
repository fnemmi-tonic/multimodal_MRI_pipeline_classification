###This is an example script that shows the use of all relevant functions###

#set the directory as your working dir: this directory must contains
#all relevant scripts
setwd("C:/my_working_dir")
#load relevant libraries and functions
library(RWeka)
library(tidyverse)
source("helper_functions.R")
source("functions_for_fitting_all_possible_combination.R")
source("fit_and_eval_all_combo.s")
#create one character vector for each modality containing the path to the directory
#where the relevant images are stored, the name of the relevant mask (that must be in the same dir)
#and a string used to search relevant images (e.g. swr for SMP processed images)
gm_vec <- c("C:/images_dir/t1", "gm_mask.nii.gz", "smoothed_normalized_gm")
fa_vec <- c("C:/images_dir/fa", "wm_mask.nii.gz", "smoothed_normalized_fa")
md_vec <- c("C:/images_dir/md", "whole_brain_mask.nii.gz", "smoothed_normalized_md")
#create the structures (a named list) that will be used for analysis
matrices <- extract_and_normalize_matrix(gm = gm_vec, fa = fa_vec, md = md_vec)
#create a dataframe with the outcome of interest (in the same order than the images)
outcome <- data_frame(outcome = c(rep("Park",26), rep("HC", 26)))
#create a vector for CV (here a 10-fold)
fold <- caret::createFolds(outcome$outcome, k = 10, list = FALSE)
#run the pipeline: here we provide the data, the label (outcome), the vector with the folding
#a vector specifing which folds we want to be treated (useful for pseudoparallelization)
#a NULL for subjects ID (see the function for information on this) and a minimum extent threshold
#for the clustering algorithm
mod_out <- fit_and_eval(matrices,outcome,fold, c(1:10), subjects_id = NULL, minimum_extent = 10)
