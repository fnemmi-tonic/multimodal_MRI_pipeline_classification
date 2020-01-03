#fit_and_eval
#this is the main function in the package. It is a wrap around the other helper function and implement 
#the pipeline from the input to the output of the final modem
#---PARAMETERS---
#list_of_modalities: list, this is the output of a call to extract_and_normalize_matrix
#outcome: character or numeric vector (transformed internally to a foator), the subjects labels
#fold_to_evaluate: nueric vector, the folding structure (e.g. the output of a call to caret::createFolds(y = outcome, k = 10, list = FALSE)) of length equal to th length of outcome
#fold_range: numeric vector, which folds should be evaluated (e.g. c(1,2,3) will evaluate only fold 1,2 and 3). Default to NULL, in this case all folds are evaluated
#subjects_id: character or numeric vector with an id for each subject, length equal to the length of outcome
#unique_identifier: string, a string to prefix to the backup RData files that will be written for each fold evaluated. Default to "BU_"
#save_relieff: boolean. Should the features surviving the relieff step be saved in the output list ? Default to TRUE
#skip relieff: boolean. Should the relieff step be skipped ? Non compatible with save_relieff = FALSE.
#previous relieff: list, the output of a previously run fit_and_eval function with save_relieff set to TRUE
#...: arguments to be passed to methods such as cluster_voxels
#---OUTPUT---
#the function output a list of length equal to the number of folds containing the following subfield:
all_coordinates, 
accuracy = accuracy, 
weights = SMO_weights,
fold_subjects = fold_subjects,
relieff_survivor = relieff_survivor

the coordinate of the selected
###clusters for each fold, a dataframe with ground truth and classification for each fold, the weights
###of the SMO model for each fold,  the subjects in the test and train set for each fold
#the selected modalities for each fold
#NB at variance with fit_and_eval_all_combo, this function does not test all possible modality
###combination, but limits itself to the subset selection algorithm

library(foreach)
fit_and_eval <- function(list_of_modalities, 
                         outcome, 
                         fold_to_evaluate, 
                         fold_range = NULL, 
                         subjects_id = NULL, 
                         unique_identifier = "BU_",
                         save_relieff = TRUE,
                         skip_relieff = FALSE,
                         previous_relieff = NULL,
                         ...) {
  
  if (length(fold_range) == 0) {up_to_fold <- 1:max(fold_to_evaluate)} else {up_to_fold <- fold_range}
  
  SMO_classifier <- make_Weka_classifier("weka/classifiers/functions/SMO")
  
  out <- foreach(fold_index = up_to_fold, .inorder = FALSE, 
                 .packages = c("tidyverse","dplyr", "CORElearn", "spatstat", "numDeriv", "quantmod", "Biocomb", "RWeka"),
                 .export = c("sd_thresholding_for_categorical_outcome_variables_vec", "select_features_relieff_derivatives_threshold_CORElearn",
                             "extract_weights_from_SMO", "SMO_classifier", "list_of_modalities", "outcome", "fold_to_evaluate")) %do% {
                               
                               all_mods_train <- list()
                               all_relieff_features <- list()
                               all_coordinates <- list()
                               
                               if (save_relieff) {relieff_survivor <- vector("list", length(list_of_modalities))
                               
                               names(relieff_survivor) <- names(list_of_modalities)} else {relieff_survivor = NULL}
                                 
                               
                               print(paste("working on TRAINING SET fold",fold_index, sep = " "))
                               for (mod in 1:length(list_of_modalities)) {
                               
                                
                                train <- list_of_modalities[[mod]]$matrix[fold_to_evaluate != fold_index,]
                                outcome_train <- outcome[fold_to_evaluate != fold_index]
                                img_dim <- list_of_modalities[[mod]]$img_dim
                                
                                name_of_mod <- names(list_of_modalities)[[mod]]
                                
                                features_type <- list_of_modalities[[mod]]$features_type
								                if (length(subjects_id) == 0) {training_subjects = NULL
                                  test_subjects = NULL} else {training_subjects = subjects_id[fold_to_evaluate != fold_index]
                                  test_subjects = subjects_id[fold_to_evaluate == fold_index]}
                                
                                if (!skip_relieff){
                                  
                                  print(paste("relieff, modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                  relieff <- threshold_variance_and_relieff (train,
                                                                             outcome_train,
                                                                             var_threshold = .25,
                                                                             outcome_variable = "outcome",
                                                                             estimator = "ReliefFequalK",
                                                                             features_type = features_type,
                                                                             thr_nps = "box_upper")
                                relieff_survivor[[mod]] <- relieff
                                
                                } else {print("Loading relieff survivors from previous analysis")
                                  relieff <- previous_relieff[[fold_index]]$relieff_survivor[[name_of_mod]]}
                                
                                if (features_type == "image"){
                                #coordinates finding 
                                  
                                  print(paste("coordinates finding, modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                  coordinates_from_features_colnames <- relieff[,-1] %>%
                                    colnames(.) %>%
                                    str_split(., pattern = "X") %>%
                                    map_chr(~`[`(.,2)) %>%
                                    as.numeric(.) %>%
                                    arrayInd(., img_dim) %>%
                                    as_data_frame(.)
                               
                                  coordinates_from_features_colnames$index <- relieff[,-1] %>%
                                    colnames(.) %>%
                                    str_split(., pattern = "X") %>%
                                    map_chr(~`[`(.,2)) %>%
                                    as.numeric(.)
                               
                                  print(paste("clustering, modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                  coordinates_from_features_colnames <- cluster_voxels(coordinates_from_features_colnames, ...)
                                
                                if(nrow(coordinates_from_features_colnames) == 0) {
                                  print("WARNING: this modality has no usable features")
                                  all_mods_train[[mod]] <- NA
                                  names(all_mods_train)[mod] <- name_of_mod 
                                  all_relieff_features[[mod]] <- NA
                                  names(all_relieff_features)[mod] <- name_of_mod 
                                  all_coordinates[[mod]] <- NA
                                  names(all_coordinates)[mod] <- name_of_mod 
                                  next
                                }
                               
                               
                                #averaging of clusters
                                df_clusters <- relieff %>%
                                 select(., -outcome) %>%
                                 mutate(subject = row_number()) %>%
                                 gather(., feature, value, -subject) %>%
                                 mutate(., index = str_split(feature, pattern = "X") %>% map(~`[`(.,2)) %>% unlist() %>% as.numeric()) %>%
                                 inner_join(., select(coordinates_from_features_colnames, index, cluster_id), by = "index") %>%
                                 group_by(subject, cluster_id) %>%
                                 summarise(mean_of_cluster = mean(value)) %>%
                                 spread(cluster_id, mean_of_cluster) %>%
                                 ungroup(.) %>%
                                 select(-subject) } else if (features_type == "nps"){
                                   df_clusters <- relieff %>% 
                                     select(-dummy)
                                   coordinates_from_features_colnames <- NA}
                                
                                relieff_features <- colnames(relieff)[-1]
                                rm(relieff)
                               
                                #change names of features for easiness of recognition later on
                                colnames(df_clusters) <- paste(name_of_mod, colnames(df_clusters),sep ="_")
                                all_mods_train[[mod]] <- df_clusters
                                names(all_mods_train)[mod] <- name_of_mod 
                                all_relieff_features[[mod]] <- relieff_features
                                names(all_relieff_features)[mod] <- name_of_mod 
                                all_coordinates[[mod]] <- coordinates_from_features_colnames
                                names(all_coordinates)[mod] <- name_of_mod }
                               
                               if(sum(is.na(all_mods_train)) != 0) { 
                                 all_mods_train <- all_mods_train[-which(is.na(all_mods_train))]}
                               
                                merged_modalities_df <- Reduce(bind_cols, all_mods_train)
                               
                                merged_modalities_df$outcome <- outcome_train
                               
                                merged_modalities_df_selected <- merged_modalities_df %>%
                                 select(., select.cfs(merged_modalities_df)$Index, outcome)
                               
                                rm(merged_modalities_df)
                               
                               #cluster and select test set
                               
                              all_mods_test <- list()
                              print("working on TEST SET")
                              for (mod in 1:length(list_of_modalities)) {
                                test <- list_of_modalities[[mod]]$matrix[fold_to_evaluate == fold_index,]
                                img_dim <- list_of_modalities[[mod]]$img_dim
                                name_of_mod <- names(list_of_modalities)[mod]
                                outcome_test <- outcome[fold_to_evaluate == fold_index]
                                features_type <- list_of_modalities[[mod]]$features_type
                                print(paste("modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                
                                if (is.na(all_relieff_features[mod])){ 
                                  print ("WARNING: this modality has no usable features")
                                  all_mods_test[[mod]] <- NA
                                  names(all_mods_test)[mod] <- name_of_mod
                                  next}
                                
                               
                               if (features_type == "image") {
                                 test_selected <- test %>%
                                 select(., all_relieff_features[[mod]]) %>%
                                 mutate(subject = row_number()) %>%
                                 gather(., feature, value, -subject) %>%
                                 mutate(., index = str_split(feature, pattern = "X") %>% map(~`[`(.,2)) %>% unlist() %>% as.numeric()) %>%
                                 inner_join(., select(all_coordinates[[mod]], index, cluster_id), by = "index") %>%
                                 group_by(subject, cluster_id) %>%
                                 summarise(mean_of_cluster = mean(value)) %>%
                                 spread(cluster_id, mean_of_cluster) %>%
                                 ungroup(.) %>%
                                 select(-subject) } else if (features_type == "nps") {
                                   test_selected <- test %>%
                                     select(., all_relieff_features[[mod]])
                                 }
                                 
                               colnames(test_selected) <- paste(name_of_mod, colnames(test_selected),sep ="_")
                               all_mods_test[[mod]] <- test_selected
                               names(all_mods_test)[mod] <- name_of_mod}
                              
                              if(sum(is.na(all_mods_test)) != 0) { 
                                all_mods_test <- all_mods_test[-which(is.na(all_mods_test))]}
                              
                              
                               merged_modalities_df_test <- Reduce(bind_cols, all_mods_test) %>%
                                 select(., head(colnames(merged_modalities_df_selected),-1))
                               
                               
                               
                               model_SMO <- SMO_classifier(as.factor(outcome) ~ ., data = merged_modalities_df_selected)
                               
                               SMO_weights <- extract_weights_from_SMO(model_SMO)
                               
                               classification <- predict(model_SMO, merged_modalities_df_test)
                               
                               
                               accuracy <- data_frame(classification = classification, ground = outcome_test)
							   fold_subjects <- list(test_subjects = test_subjects, training_subjects = training_subjects)
                               out <- list(all_coordinates, 
                                           accuracy = accuracy, 
                                           weights = SMO_weights,
                                           fold_subjects = fold_subjects,
                                           relieff_survivor = relieff_survivor)
                               backup_file_name <- paste(unique_identifier, "_fold_", fold_index,".RData", sep = "")
                               save(file = backup_file_name, out)
                               out
                               
                             }
  
}
  
  