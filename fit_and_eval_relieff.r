#this function is a wrap around most of the other functions
#it takes as input the list output of the function extract_and_normalize_matrix.R,
#the label of the outcome in the form of a vector
#a vector of folding label (a number from 1 to N folds for each subject)
#a vector stating which fold to evaluate (useful for pseudoparallelization)
#a char vector of the same size of outcome with a label for each subject (useful to keep track of which
###subject are wrongly classified) or a NULL
#a unique identifier that is used to save a backup copy of the results of each fold
#the function output a list of length N of fold that contains the coordinate of the selected
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
  print("I have at least entered this fucking function")
  if (length(fold_range) == 0) {up_to_fold <- 1:max(fold_to_evaluate)} else {up_to_fold <- fold_range}
  print("line 29 so far so good")
  SMO_classifier <- make_Weka_classifier("weka/classifiers/functions/SMO")
  print("line 31 so far so good")
  out <- foreach(fold_index = up_to_fold, .inorder = FALSE, 
                 .packages = c("tidyverse","dplyr", "CORElearn", "spatstat", "numDeriv", "quantmod", "Biocomb", "RWeka"),
                 .export = c("sd_thresholding_for_categorical_outcome_variables_vec", "select_features_relieff_derivatives_threshold_CORElearn",
                             "extract_weights_from_SMO", "SMO_classifier", "list_of_modalities", "outcome", "fold_to_evaluate")) %do% {
                               print("I have entered the foreach loop")
                               all_mods_train <- list()
                               all_relieff_features <- list()
                               all_coordinates <- list()
                               
                               if (save_relieff) {relieff_survivor <- vector("list", length(list_of_modalities))
                               
                               names(relieff_survivor) <- names(list_of_modalities)} else {relieff_survivor = NULL}
                                 
                               
                               print(paste("working on TRAINING SET fold",fold_index, sep = " "))
                               for (mod in 1:length(list_of_modalities)) {
                               print("I have entered the mod loop") 
                                
                                train <- list_of_modalities[[mod]]$matrix[fold_to_evaluate != fold_index,]
                                outcome_train <- outcome[fold_to_evaluate != fold_index]
                                img_dim <- list_of_modalities[[mod]]$img_dim
                                print("Imma bout to read the name_of_mod")
                                name_of_mod <- names(list_of_modalities)[[mod]]
                                print(str_glue("I fucking read the name of mod, that by the way is {name_of_mod}"))
                                features_type <- list_of_modalities[[mod]]$features_type
								                if (length(subjects_id) == 0) {training_subjects = NULL
                                  test_subjects = NULL} else {training_subjects = subjects_id[fold_to_evaluate != fold_index]
                                  test_subjects = subjects_id[fold_to_evaluate == fold_index]}
                                
                                if (!skip_relieff){
                                  print("started the fucking relieff")
                                  print(paste("relieff, modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                  relieff <- threshold_variance_and_relieff (train,
                                                                             outcome_train,
                                                                             var_threshold = .25,
                                                                             outcome_variable = "outcome",
                                                                             estimator = "ReliefFequalK",
                                                                             features_type = features_type,
                                                                             thr_nps = "box_upper")
                                relieff_survivor[[mod]] <- relieff
                                print("ended relieff")
                                } else {print("Loading relieff survivors from previous analysis")
                                  relieff <- previous_relieff[[fold_index]]$relieff_survivor[[name_of_mod]]}
                                
                                if (features_type == "image"){
                                #coordinates finding 
                                  print("entered image modality branch")
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
                                 select(-subject) } else if (features_type == "nps"){print("entered nps modality branch")
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
  
  