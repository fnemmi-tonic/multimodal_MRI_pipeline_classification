# multimodal_MRI_pipeline_classification
scripts and helper functions of a multimodal MRI pipeline for classification (binary and multiclass problems)

-helper_functions.r contains all companion functions needed to run the pipeline except for the functions needed to extensively search the modalities space for the best modalities

-functions_for_fitting_all_possible_combination.r contains the functions needed to extensively search the modalities space for the best modalities

-fit_and_eval_all_combo.r is the main function that wrap the all pipeline

-evaluation_of_model_and_cluster_writing.r contains function to evaluate the model performance and to write the most discriminative clusters