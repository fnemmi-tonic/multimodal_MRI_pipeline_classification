#load library
library(stringr)
library(tidyverse)
library(purrr)
library(oro.nifti)
library(neurobase)
library(CORElearn)
library(numDeriv)
library(quantmod)


#residuals_on_a_vector
#returns the residuals of a linear model fitted using nuisance variables as independent variable
##and the variable to be denoised as dependent variable
#---PARAMETERS---
#vec: numeric vector, the variable to be denoised
#iv_as_a_dataframe: nuisance variables as data_frame or data.frame
#formula: a formula to fit a linear model. if NULL the formula is automatically set as "~.". Default to null
#---OUTPUT---
##res : numeric vector with the residuals
residuals_on_a_vector <- function(vec, iv_as_a_dataframe, formula = NULL) {
  if (is.null(formula)) {formula <- as.formula("dependent ~ .")} else {formula <- as.formula(paste("dependent ~ ",formula))}
  df <- iv_as_a_dataframe
  df$dependent <- vec
  model <- lm(formula = formula, data = df)
  res <- residuals(model)
  return(res)
}

#residuals_on_a_dataframe
#adaptation of residuals_on_a_vector for whole dataframe
#---PARAMETERS---
#df: a data_frame or data.frame with rows = number of subjects and columns = variables to be denoised
#iv_as_a_dataframe: nuisance variables as data_frame or data.frame
#formula: a formula to fit a linear model. if NULL the formula is automatically set as "~.". Default to null
#---OUTPUT---
##out : denoised data_frame
residuals_on_a_dataframe <- function(df, iv_as_a_dataframe, formula = NULL) {
  
  out <- map(df,residuals_on_a_vector, iv_as_a_dataframe, formula)
  return(as_data_frame(out))
  
}

#extract_weights_from_SMO
#extract weights and variable names from a fitted SMO model (discrimination)
#---PARAMETERS---
#model: a fitted SMO model
#---OUTPUT---
##df: a data_frame with columns "names" and "weights' containing the variables names and weights
extract_weights_from_SMO <- function(model) {
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  
  raw_output <- capture.output(model)
  trimmed_output <- raw_output[-c(1:11,(length(raw_output) - 4): length(raw_output))]
  df <- data_frame(features_name = vector(length = length(trimmed_output) + 1, "character"), 
                   features_weight = vector(length = length(trimmed_output) + 1, "numeric"))
  
  if (length(trimmed_output) == 0) {
  df$features_name <- NA
  df$features_weight <- NA} else {
  
  for (line in 1:length(trimmed_output)) {
    
    
    string_as_vector <- trimmed_output[line] %>%
      str_split(string = ., pattern = " ") %>%
      unlist(.)
    
    
    numeric_element <- trimmed_output[line] %>%
      str_split(string = ., pattern = " ") %>%
      unlist(.) %>%
      as.numeric(.)
    
    position_mul <- string_as_vector[is.na(numeric_element)] %>%
      str_detect(string = ., pattern = "[*]") %>%
      which(.)
    
    numeric_element <- numeric_element %>%
      `[`(., c(1:position_mul))
    
    text_element <- string_as_vector[is.na(numeric_element)]
    
    
    there_is_plus <- string_as_vector[is.na(numeric_element)] %>%
      str_detect(string = ., pattern = "[+]") %>%
      sum(.)
    
    if (there_is_plus) { sign_is <- "+"} else { sign_is <- "-"}
    
    
    
    feature_weight <- numeric_element[!is.na(numeric_element)]
    
    if (sign_is == "-") {df[line, "features_weight"] <- feature_weight * -1} else {df[line, "features_weight"] <- numeric_element[!(is.na(numeric_element))]}
    
    df[line, "features_name"] <- paste(text_element[(position_mul + 1): length(text_element)], collapse = " ")
    
  }
  
  intercept_line <- raw_output[length(raw_output) - 4]
  
  
  there_is_plus_intercept <- intercept_line %>%
    str_detect(string = ., pattern = "[+]") %>%
    sum(.)
  
  if (there_is_plus_intercept) { intercept_sign_is <- "+"} else { intercept_sign_is <- "-"}
  
  numeric_intercept <- intercept_line %>%
    str_split(string = ., pattern = " ") %>%
    unlist(.) %>%
    as.numeric(.) %>%
    `[`(., length(.))
  
  df[nrow(df), "features_name"] <- "intercept"
  
  if (intercept_sign_is == "-") {df[nrow(df), "features_weight"] <- numeric_intercept * -1} else {df[nrow(df), "features_weight"] <- numeric_intercept}
  
  options(warn = oldw)
  
  
  df <- df %>%
    arrange(desc(abs(features_weight)))}
  
  return(df)
}

#extract_weights_from_SMOreg
#extract weights and variable names from a fitted SMO model (regression)
#---PARAMETERS---
#model: a fitted SMO model
#---OUTPUT---
##df: a data_frame with columns "names" and "weights' containing the variables names and weights
extract_weights_from_SMOreg <- function(model) {
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  
  raw_output <- capture.output(model)
  trimmed_output <- raw_output[-c(1:3,(length(raw_output) - 4): length(raw_output))]
  df <- data_frame(features_name = vector(length = length(trimmed_output) + 1, "character"), 
                   features_weight = vector(length = length(trimmed_output) + 1, "numeric"))
  
  for (line in 1:length(trimmed_output)) {
    
    
    string_as_vector <- trimmed_output[line] %>%
      str_split(string = ., pattern = " ") %>%
      unlist(.)
    
    
    numeric_element <- trimmed_output[line] %>%
      str_split(string = ., pattern = " ") %>%
      unlist(.) %>%
      as.numeric(.)
    
    position_mul <- string_as_vector[is.na(numeric_element)] %>%
      str_detect(string = ., pattern = "[*]") %>%
      which(.)
    
    numeric_element <- numeric_element %>%
      `[`(., c(1:position_mul))
    
    text_element <- string_as_vector[is.na(numeric_element)]
    
    
    there_is_plus <- string_as_vector[is.na(numeric_element)] %>%
      str_detect(string = ., pattern = "[+]") %>%
      sum(.)
    
    if (there_is_plus) { sign_is <- "+"} else { sign_is <- "-"}
    
    
    
    feature_weight <- numeric_element[!is.na(numeric_element)]
    
    if (sign_is == "-") {df[line, "features_weight"] <- feature_weight * -1} else {df[line, "features_weight"] <- numeric_element[!(is.na(numeric_element))]}
    
    df[line, "features_name"] <- paste(text_element[(position_mul + 1): length(text_element)], collapse = " ")
    
  }
  
  intercept_line <- raw_output[length(raw_output) - 4]
  
  
  there_is_plus_intercept <- intercept_line %>%
    str_detect(string = ., pattern = "[+]") %>%
    sum(.)
  
  if (there_is_plus_intercept) { intercept_sign_is <- "+"} else { intercept_sign_is <- "-"}
  
  numeric_intercept <- intercept_line %>%
    str_split(string = ., pattern = " ") %>%
    unlist(.) %>%
    as.numeric(.) %>%
    `[`(., length(.))
  
  df[nrow(df), "features_name"] <- "intercept"
  
  if (intercept_sign_is == "-") {df[nrow(df), "features_weight"] <- numeric_intercept * -1} else {df[nrow(df), "features_weight"] <- numeric_intercept}
  
  options(warn = oldw)
  
  
  df <- df %>%
    arrange(desc(abs(features_weight)))
  
  return(df)
}

#normalize_matrix_range
#normalize matrix by rescaling from 0 and 1
#---PARAMETERS---
#matrix: the matrix to be scaled
#---OUTPUT---
##new_mat: the range normalized matrix
#normalize matrix range
normalize_matrix_range <- function(matrix) {
  range_mat <- range(matrix)
  new_mat <- (matrix - range_mat[1])/(range_mat[2] - range_mat[1])
  return(new_mat)
}


normalize_vector_range <- function(vec){
  range_vec <- range(vec)
  new_vec <- (vec - range_vec[1])/(range_vec[2] - range_vec[1])
  return(new_vec)
}


normalize_matrix_range_by_columns <- function(mat){
  new_mat <- apply(mat, 2, normalize_vector_range)
  return(new_mat)
}


#reshape_images_for_pipeline
#Reshape two or more 3D images in nifti or nifti gz format to a matrix with number of rows = number of images
## and number of columns = number of voxels in the images. If the features are 2D (csv file), read the csv file and create the 
##expected matrix
#---PARAMETERS---
#image_dir: the directory where all the images (of the same dimension) are stored
#mask : string, the name of the image to use as mask (it should be in the same directory as image_dir). 
##If the features are 2D (features_type = "nps") enter NULL. default to NULL
#pattern_for_search: string, the pattern to search for listing all the images (e.g. if you have smoothed normalized
##images from SPM this could be "sw") or the csv file if  
#features_type: string, the type of features entered (accepted are "nps" and "image")
#delim: only useful for csv file (2D features)
#subjects_index: indexes of subject to be included, if NULL all subjects are included. Default to NULL
#---OUTPUT---
##new_mat: the range normalized matrix
reshape_images_for_pipeline <- function (image_dir, 
                                         mask = NA, 
                                         pattern_for_search,
                                         features_type,
                                         delim = ";",
                                         subjects_index = NULL) {
  
  this_wd <- getwd()
  setwd(image_dir)
  
  if (features_type == "image"){
  list <- dir(pattern = pattern_for_search)
  number_of_subjects <- length(list)
  
  if (is.null(subjects_index)) {subjects <- 1:number_of_subjects} else {subjects <- subjects_index}
  
  number_of_subjects <- length(subjects)
  
  
  
  
  print("reading mask")
  mask_struct <- readNIfTI2(mask)
  mask_sparse <- which(mask_struct@.Data > 0)
  dimension <- length(mask_sparse)
  n_by_v_matrix <- matrix(data = NA, nr = length(subjects), nc = Reduce(`*`, dimension))
  
  for (image_index in 1:length(subjects)) {
    print(paste("reading and processing subject", subjects_index[image_index], "file is", list[subjects[image_index]]))
    img_struct <- readNIfTI2(list[subjects[image_index]])
    n_by_v_matrix[image_index,] <- img_struct@.Data[mask_sparse]
    
  }
  
  colnames(n_by_v_matrix) <- paste("X", mask_sparse, sep = "")
  setwd(this_wd)
  n_by_v_matrix <- as_data_frame(n_by_v_matrix)
  list(n_by_v_matrix = n_by_v_matrix, 
       dim_img = dim(mask_struct@.Data), 
       img_struct = mask_struct,
        features_type = features_type)}
  
  else if (features_type == "nps") {
    list <- dir(pattern = pattern_for_search)
    n_by_v_matrix <- read_delim(list[1], delim = delim) %>% 
      as.matrix()
    if (is.null(subjects_index)) {subjects <- 1:nrow(n_by_v_matrix)} else {subjects <- subjects_index}
    n_by_v_matrix <- n_by_v_matrix[subjects,]
    setwd(this_wd)
    list(n_by_v_matrix = as_data_frame(n_by_v_matrix),
         dim_img = dim(n_by_v_matrix),
         img_struct = NULL,
         features_type = features_type)}
    
    
    
  }
  


#extract_and_normalize_matrix
#Wrapper around reshape_images_for_pipeline
#---PARAMETERS---
#...: a number of named vectors with names equal to the parameters of reshape_images_for_pipeline
#subjects_index: indexes of subject to be included, if NULL all subjects are included. Default to NULL
#---OUTPUT---
##a list of matrices of length equal to the input vector
extract_and_normalize_matrix <- function(...,
                                         subjects_index = NULL) {
  
  arguments <- list(...)
  
  for (arg in 1:length(arguments)) {
    
    info <- reshape_images_for_pipeline(arguments[[arg]][1], 
                                        arguments[[arg]][2], 
                                        arguments[[arg]][3],
                                        arguments[[arg]][4],
                                        arguments[[arg]][5],
                                        subjects_index)
    matrix <- info$n_by_v_matrix
    img_dim <- info$dim_img
    feat_type <- info$features_type
    out <- list(matrix = matrix, img_dim = img_dim, features_type = feat_type)
    assign(names(arguments)[arg], out)
  }
  
  return(mget(names(arguments)))
}



#sd thresholding for categorical outcome variable
library(magrittr)
library(dplyr)

#var_vectorized
#calculate variance of a matrix by columns in a vectorized way
#---PARAMETERS---
#mat: the numeric matrix
#---OUTPUT---
#a vector with variance for each columns
var_vectorized <- function(mat) {
  
  mean_r <- colSums(mat)/nrow(mat)
  mat_sub_mean <- (t(mat) - mean_r)^2
  var_r <- rowSums(mat_sub_mean)/(ncol(mat_sub_mean)-1)
}

#sd_thresholding_for_categorical_outcome_variables_vec
#retains only those features that have variance higher than a certain threshold
#---PARAMETERS---
#df: a data_frame of features
#quant: numeric, the threshold
#---OUTPUT---
#a data_frame wich include only the retained variables
sd_thresholding_for_categorical_outcome_variables_vec <- function(df, quant) {
  
  
  var_each_features <- var_vectorized(df)
  
  var_thr <- quantile(var_each_features, quant)
  
  features_above_treshold <- which(var_each_features > var_thr)
  
  output_df <- df %>%
    select(., features_above_treshold)
  
  return(output_df)
}



#relieff feat selection
#calculate the derivative of a GRAPHICAL function (i.e. a function for which you have x and y)
#---PARAMETERS---
#x: numeric vector, coordinates on the x axis
#y: numeric vector, coordinates on the y axis
#---OUTPUT---
#the derivative
deriv <- function(x,y) {
  output <- diff(y)/diff(x)
  return(output)
}

#middle pts
#calculate the middle points of lagged difference of a numerical series
#---PARAMETERS---
#x: the numerical series whose middle points should be calculated
#---OUTPUT---
#a numerical vector with the middle points
middle_pts <- function(x){
  pts <- x[-1] - diff(x) / 2
  return(pts)
}

#calculate_features_threshold_based_on_second_derivative
#calculate the reliefF treshold based on the scree plot methods (i.e. find the minimum of the second derivatives
##of the orderd reliefF weights)
#---PARAMETERS---
#x: numeric vector, the ordinal position of the weight (i.e. seq(1,length(y)))
#y: ordered relieff weights
#to_plot: boolean, should the weigths and the threshold be plotted ? default to TRUE
#---OUTPUT---
#a named list containing the smoothed second derivative, the original second derivative and the threshold
calculate_features_threshold_based_on_second_derivative <- function(x,y, to_plot = TRUE) {
  smoothed_y <- loess(y ~ x,
                      data.frame(y = y, x = x, model = T), span = .1) %>% predict(.)
  
  second_d <- deriv(middle_pts(x), deriv(x, smoothed_y))
  smooth_second_d <- loess(second_d ~ midpts,
                           data.frame(second_d = second_d, midpts = middle_pts(middle_pts(x))), model = T, span = .1)
  otp <- predict(smooth_second_d)
  thr <- y[findValleys(otp)[1]]
  
  if(to_plot) {
    plot(x,y)
    points(findValleys(otp)[1],y[findValleys(otp)[1]], pch = 15, cex = 2, col = "magenta")
  }
  return(list(smoothed_second_derivatives = smooth_second_d, second_derivative_values = otp, threshold = thr))
}

#calculate_features_threshold_based_on_second_derivative
#calculate the reliefF weights for the features in a dateframe and return a dataframe containing only the features
##with reliefF weight higher than the calculated threshold (see function above)
#---PARAMETERS---
#df: data.frame or (preferred) data_frame structures with nrow = number of subjects and ncol = number of features + 1
##the last column of the dataframe should contain the outcome
#outcome: string, the name of the last column
#estimator: one of the accepted estimator of CORElearn::attrEval
#---OUTPUT---
#a data_frame containing only the features surviving the reliefF threshold plus the outcome
select_features_relieff_derivatives_threshold_CORElearn <- function(df, outcome, estimator) {
  
  
  
  print("Performing relieff algorithm")
  
  rf_weights <- attrEval(formula = ncol(df), df, estimator = estimator)
  
  print("Done relieff algorithm - calculating threshold")
  
  ordered_weights <- sort(rf_weights, decreasing = TRUE)
  
  thr <- calculate_features_threshold_based_on_second_derivative(seq(1,length(ordered_weights)), ordered_weights)$threshold
  
  selected_weights <- names(rf_weights) [rf_weights > thr]
  
  output <- df %>%
    select(one_of(c(outcome,selected_weights)))
  
}





cluster_voxels <- function(coordinates_table, minimum_extent = 10, distances = c(1.8,3,4), n_passes = 3){
  
  if (length(distances) != n_passes) {
    warning("you have not specified the correct number of distances")
  }
  
  condition = TRUE
  index = 0
  
  while (condition == TRUE && index < n_passes) {
    
    index = index + 1
    print(index)
    if (!is.numeric(range(coordinates_table[,1])) || length(range(coordinates_table[,1])) < 2 || diff(range(coordinates_table[,1])) == 0 ||
        !is.numeric(range(coordinates_table[,2])) || length(range(coordinates_table[,2])) < 2 || diff(range(coordinates_table[,2])) == 0 ||
        !is.numeric(range(coordinates_table[,3])) || length(range(coordinates_table[,3])) < 2 || diff(range(coordinates_table[,3])) == 0) {break}
    
    bb <- box3(range(coordinates_table[,1]),
               range(coordinates_table[,2]), range(coordinates_table[,3]))
    object.pp3 <- pp3(coordinates_table$V1,
                      coordinates_table$V2, coordinates_table$V3, bb)
    object.pp3_labelled <- connected(object.pp3, R = distances[index])
    coordinates_table$cluster_id <- marks(object.pp3_labelled)
    cluster_extent <- table(coordinates_table$cluster_id)
    cluster_extent <- data_frame(extent = as.vector(cluster_extent), cluster_id = names(cluster_extent))
    cluster_extent <- left_join(coordinates_table, cluster_extent, by = "cluster_id")
    pass_name <- paste("cluster_to_retain_at_pass_", index, sep ="")
    assign(pass_name, cluster_extent %>% filter(extent >= minimum_extent))
    if (sum(table(coordinates_table$cluster_id) < minimum_extent) == 0) {condition = FALSE}
    coordinates_table <- cluster_extent %>% 
      filter(extent < minimum_extent) %>% 
      select(V1, V2, V3, index)
  }
  
  all_passes <- mget(ls(patt = "cluster_to_retain"))
  
  if(length(all_passes) == 1) {return(all_passes[[1]])}
  
  all_passes <- map2(all_passes,seq(1,length(all_passes)), ~ mutate(.x,set = .y)) %>%
    Reduce(bind_rows, .) %>%
    mutate(tt = paste(cluster_id, set, sep = "_")) %>% 
    arrange(tt) %>% 
    mutate(old_ind = as.numeric(as.factor(tt))) %>%
    select(-extent, -set,-tt,-cluster_id, cluster_id = old_ind)
  
  return(all_passes)
}


create_n_balanced_folds <- function(original_fold, outcome, n_folds, maxIter = 10000) {
  
  opt_list <- list()
  fisher_ps <- vector()
  iteration_counter <- 0
  index = 0
  test_folds <- n_folds
  while (test_folds != 0) {
    this_fold <- sample(fold)
    fisher_p <- fisher.test(table(this_fold, outcome))$p.value
    eq_test <- sum(unlist(map2(list(this_fold), opt_list, identical)))
    if (fisher_p > .2 && eq_test == 0) {
      index = index + 1
      opt_list[[index]] <- this_fold
      test_folds <- test_folds - 1
      fisher_ps[index] <- fisher_p
    }
    iteration_counter <- iteration_counter + 1
    
    if (iteration_counter == maxIter) {
      warning("maxIter reached")
      break
    }
  }
  
  return(list(folds = opt_list, fisher = fisher_ps))
  
}



threshold_variance_and_relieff <- function(train_features, 
                                           train_outcome,
                                           var_threshold = .25,
                                           outcome_variable = "outcome",
                                           estimator, 
                                           features_type,
                                           thr_nps = "box_upper"){
  if (estimator == "ReliefFequalK") {train_outcome <- as.factor(train_outcome)}
  if (features_type == "nps"){
    var_thresholded <- sd_thresholding_for_categorical_outcome_variables_vec(train_features, var_threshold)
    var_thresholded$outcome <- train_outcome
    df <- var_thresholded
    rf_weights <- attrEval(formula = ncol(df), df, estimator = estimator)
    switch (thr_nps,
            upper_whisker <- {relief_thr <- boxplot(rf_weights, plot = FALSE)$stats[5]},
            box_upper = {relief_thr <- boxplot(rf_weights, plot = FALSE)$stats[4]},
            median = {relief_thr <- boxplot(rf_weights, plot = FALSE)$stats[3]},
            box_lower = {relief_thr <- boxplot(rf_weights, plot = FALSE)$stats[2]},
            lower_whisker = {relief_thr <- boxplot(rf_weights, plot = FALSE)$stats[1]})
    surviving_features_name <- names(rf_weights[rf_weights >= relief_thr])
    surviving_features <- train_features %>% select(one_of(surviving_features_name)) %>%
      mutate(dummy = 1) %>%
      select(dummy,one_of(colnames(.)))} else if (features_type == "image") {
      var_thr <- sd_thresholding_for_categorical_outcome_variables_vec(train_features, var_threshold)
      var_thr$outcome <- train_outcome
      surviving_features <- select_features_relieff_derivatives_threshold_CORElearn(var_thr, outcome_variable, estimator = estimator)}
  
  return(surviving_features)
}



