library(oro.nifti)
library(neurobase)
library(tidyverse)
library(stringr)

#output_selected_clusters
#this function output, for each fold, the selected clusters for each modality as nifti images
#---PARAMETERS---
#output_list: the list that is output of a call of the fit_and_eval function
#space_defining_image : the absolute path to a nii or a nii.gz file defining the space in which to write the clusters
#output_name: string, prefix of the filenames of the nifti images to be written
output_selected_clusters <- function(output_list, space_defining_image, output_name) {
  
  img <- readNIfTI(space_defining_image)
  
  img_dim <- dim(img@.Data)
  
  folds <- length(output_list)
  
    
  for (fold_index in 1:folds) {
  
  modalities <- names(output_list[[fold_index]][[1]])
  
  modalities_number <- length(modalities)
    
    for (modality_index in 1:modalities_number) { 
      
      this_modality <- modalities[modality_index]
      
      clusters_number <- sum(str_detect(output_list[[fold_index]]$weights$features_name, this_modality))
      
      if(clusters_number == 0){
        next}
      
      suffix <- paste("fold", fold_index, this_modality, sep = "_")
      
	  
	  if (this_modality == "dorsal_attentional") {
		this_modality_clusters <- str_subset(output_list[[fold_index]]$weights$features_name, this_modality) %>%
        str_split(.,"_") %>%
        map_chr(~`[`(.,3)) %>%
        as.numeric(.)} else {
	  
      this_modality_clusters <- str_subset(output_list[[fold_index]]$weights$features_name, this_modality) %>%
        str_split(.,"_") %>%
        map_chr(~`[`(.,2)) %>%
        as.numeric(.)}
      
      selected_clusters <- `[[`(output_list[[fold_index]][[1]], this_modality) %>%
      filter(., cluster_id %in% this_modality_clusters)
    
      new_img <- img
      new_img@.Data <- array(0,img_dim)
    
      for (cl in this_modality_clusters) {
        new_img@.Data[as.matrix(selected_clusters[selected_clusters$cluster_id == cl, 1:3])] <- cl
        
      }
    
      img_name <- paste(output_name,suffix, sep = "_")
    
      writeNIfTI(new_img, paste0(img_name), verbose=TRUE)
      
    }
    
    
  }
  
}
#evaluators_calculation
#helper function to calculate several performance indexes
#---PARAMETERS---
#df: data_frame or data.frame with columns "classification" (label from the fitted model) and "ground" (ground truth)
evaluators_calculation <- function(df) {
  
  performance_table <- table(df$classification, df$ground)
  
  accuracy <- sum(diag(performance_table))/sum(performance_table)
  sensitivity <- performance_table[2,2]/sum(performance_table[,2])
  specificity <- performance_table[1,1]/sum(performance_table[,1])
  F1 <- (2*performance_table[2,2])/((2*performance_table[2,2]) + performance_table[1,2] + performance_table[2,1])
  
  evaluation <- c(accuracy = accuracy, sensitivity = sensitivity, specificity = specificity, F1 = F1)
  
  return(evaluation)
  
}

#evaluate_model_merged
#this function take as input the output of the function fit_and_eval and output
#several performance indexes for the stacked folds
#---PARAMETERS---
#output_classification: list, the output of a call to fit_and_eval
evaluate_model_merged <- function(output_classification) {
  
  all_folds <- output_classification %>%
    map(~`$`(.,"accuracy")) %>%
    Reduce(bind_rows,.)
  
  output <- evaluators_calculation(all_folds)
  
  return(output)
}

#evaluate_model_merged
#this function take as input the output of the function fit_and_eval and output
#several performance indexes averaged by fold
#---PARAMETERS---
#output_classification: list, the output of a call to fit_and_eval
evaluate_model_fold <- function(output_classification) {
  
  each_folds <- output_classification %>%
    map(~`$`(.,"accuracy"))
  
  each_eval <- lapply(each_folds, evaluators_calculation)
  
  each_accuracy <- sapply(each_eval,`[`,"accuracy")
  each_accuracy [is.nan(each_accuracy)] <- 0
  mean_accuracy <- mean(each_accuracy)
  
  each_sensitivity <- sapply(each_eval,`[`,"sensitivity")
  each_sensitivity [is.nan(each_sensitivity)] <- 0
  mean_sensitivity <- mean(each_sensitivity)
  
  each_specificity <- sapply(each_eval,`[`,"specificity")
  each_specificity [is.nan(each_specificity)] <- 0
  mean_specificity <- mean(each_specificity)
  
  each_F1 <- sapply(each_eval,`[`,"F1")
  each_F1 [is.nan(each_F1)] <- 0
  mean_F1 <- mean(each_F1)
  
  output <- c(accuracy = mean_accuracy, sensitivity = mean_sensitivity, specificity = mean_specificity, F1 = mean_F1)
  
  return(output)
}

