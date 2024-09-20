library(tidyverse)
#library(parallel)
#library(mlr)
#library(pbapply)


# file to hold the code used for predictions accessing web application

# This script is meant to take the outputs of Diamond and shape them for ML
# Choice for metrics should be included in the config files for now 
# Make sure to adjust file inputs based on working dir (Especially the config file)


reformat_diamond <- function(diamond_data, model_config) {
    # This script is meant to take the outputs of Diamond and shape them for ML
    # Choice for metrics should be included in the config files for now 
    # Make sure to adjust file inputs based on working dir (Especially the config file)

    split_conds <- function(my_dat) {
        # splits conditions (can be replaced with group_split)
        # sets up a list of tibbles containing the data for specific conditions
        empty_list <- vector(mode = "list", length = length(model_config$conditions))
        # ^ not good because it references model_config which is not an input -DYZ
        
        # for each condition, subset data
        for (i in 1:length(model_config$conditions)) {
            conds <- model_config$conditions[i]
            names(empty_list)[i] = conds
            sub_dat = my_dat %>% filter(Condition == conds)
            empty_list[[i]] = sub_dat
            All_conds = empty_list
        }
        return(All_conds)
    } # END

    split_tps <- function(tp_dat) {
    
        const_tp = c("Time point 1") # hard-code constant time point string
        # initialize empty vector for time points subsets
        tp_list = vector(mode = 'list', length = length(model_config$conditions))
        
        # names of the conditions match those of the condition subsets
        names(tp_list) <- names(tp_dat)
        
        # loop through conditions
        for (i in 1:length(model_config$conditions)){
            dat = tp_dat[[i]] # enter into the list to get the condition subset data
            
            if (names(tp_dat)[i] %in% c("Butyrate","Valerate","Cholesterol","HighCholesterol")) {
            # get the last time point unique to that subcondition
            # ^ however, what if different drugs were tested differently?
            term_tp = tail(unique(dat$TimePoint), 1)
            cond_meas = c("OD")
            new_dat = dat %>% filter(TimePoint %in% c(const_tp,term_tp) & Measurement == cond_meas)
            
            } else if (names(tp_dat)[i] %in% c("LowpH")){
            term_tp = c("Time point 4")
            cond_meas = c("OD")
            new_dat = dat %>% filter(TimePoint %in% c(const_tp,term_tp) & Measurement == cond_meas)
            
            } else if(names(tp_dat)[i] %in% c("Macrophage")){
            term_tp = c("Time point 4")
            cond_meas = c("Lum")
            new_dat = dat %>% filter(TimePoint == term_tp & Measurement == cond_meas)
            
            } else if(names(tp_dat)[i] %in% c("Nitrate")){
            term_tp = c("Time point 3")
            cond_meas = c("Lum")
            new_dat = dat %>% filter(TimePoint %in% c(const_tp,term_tp) & Measurement == cond_meas)
            
            } else if (names(tp_dat)[i] %in% c("RichGrowth")){
            #term_tp = tail(unique(dat$TimePoint), 1)
            # ^ should be hardcoded to Time point 3, to ensure consistency between data
            term_tp <- "Time point 3"
            cond_meas = c("OD")
            new_dat = dat %>% filter(TimePoint == term_tp & Measurement == cond_meas)
            
            } else if (names(tp_dat)[i] %in% c("HighCholesterolLowpH")){
            term_tp = tail(unique(dat$TimePoint), 1)
            cond_meas = c("Lum")
            new_dat = dat %>% filter(TimePoint %in% c(const_tp,term_tp) & Measurement == cond_meas)
            }
            tp_list[[i]] = new_dat
        } 
        return(tp_list)
    } # END

    #adjust the names of each chosen conditions based common names (Can add to this if more conditions come up)
    change_Condname <- function(new_tp){
        for (i in 1:length(new_tp)){
            if (names(new_tp)[i] == "Cholesterol"){
            new_tp[[i]]$Condition = "lowcholesterol"
            } else if (names(new_tp)[i] == "LowpH"){
            new_tp[[i]]$Condition = "acidic"
            } else if (names(new_tp)[i] == "Macrophage"){
            new_tp[[i]]$Condition = "intracellular"
            } else if (names(new_tp)[i] == "Nitrate"){
            new_tp[[i]]$Condition = "dormancy"
            } else if (names(new_tp)[i] == "RichGrowth"){
            new_tp[[i]]$Condition = "standard"
            } else if(names(new_tp)[i] == "HighCholesterol"){
            new_tp[[i]]$Condition = "cholesterol"
            } else if(names(new_tp)[i] == "HighCholesterolLowpH"){
            new_tp[[i]]$Condition = "acidichighcholesterol"
            }
        }
        
        return(new_tp)
    } # END


    # adjusts for timepoint name (name it constant or terminal tp)
    change_TPname <- function(last_dat){
    
    tp_ = function(tp_name){
        if ( tp_name == "Time point 1"){
        curr_tp  = "Constant"
        } else {
        curr_tp = "Terminal"
        }
        curr_tp
    }
    for (i in 1:length(last_dat)){
        if (last_dat[[i]]$Condition[1] %in% c("intracellular","standard")){
        last_dat[[i]]$TimePoint = "ConstantTerminal"
        } else {
        last_dat[[i]]$TimePoint = map_chr(last_dat[[i]]$TimePoint, tp_)
        }
    }
    last_dat
    }

    # Obtain number of replicates
    # extracting the number of sample that went into each average
    Included_reps <- diamond_data$Included

    # Pseudo way of finding number of samples (R stands for Rep)
    num_samples <- sapply(Included_reps, str_count, "R")

    # Averaged samples represent the number of samples that went to calculating an average value
    # for each metric
    # add new column for num_samples
    diamond_data <- diamond_data %>% mutate(Averaged_samples = num_samples) 
    all_feats = names(diamond_data) # get feature names

    # generate names of all features to extract (match case sensitive)
    combined_feats <- c(model_config$meta_feats, model_config$features)
    stdev_feats <- paste("Stdev", model_config$features, sep = "")
    stat_feats <- c(stdev_feats, "Averaged_samples")
    final_feats <- c(combined_feats, stat_feats)


    # PROCESS DATA =================================================================
    subset_data <- diamond_data %>% 
        select(all_of(final_feats)) %>% 
        filter(Condition %in% model_config$conditions)

    # split (subset) data by condition and put into a list named with condition
    chosen_conds <- split_conds(subset_data) # custom function

    # can modify these time points based on more conditions
    chosen_tps = split_tps(chosen_conds)
    # ^ time points are hardcoded in. This may be the source of the error

    name_cond = change_Condname(chosen_tps)
    name_tp = change_TPname(name_cond)

    concat_alldata = bind_rows(name_tp) %>% select(-c(Measurement))
    final_cols = names(concat_alldata)[5:ncol(concat_alldata)]
    wide_format = concat_alldata %>% pivot_wider(names_from = c(Condition,TimePoint),
                                                values_from = (all_of(final_cols)))

    return(wide_format)
}

get_singles_from_combos <- function(combos_list) {
  # isolates a vector of single drugs from a vector or list of combos
  singles_list <- c() # initialize empty vector 
  for (combo in combos_list) {
    singles <- c(unlist(str_split(combo, "\\+")))
    singles_list <- append(singles_list, singles)
  }

  return(unique(sort(singles_list)))
}


generate_combos <- function(drug_list, nways, drug_pairs_exclude, labeled_data) {
  # inputs:
  #   drug_list: vector of drug abbreviations e.g. c("BDQ", "CLZ", "DEL")
  #   nways: vector of n-degree of higher orders to generate e.g. c(3, 4)
  #   drug_pairs_exclude: list of named vectors containing drugs to not pair

  combo_contains_excluded_pair <- function(combo, pairs_exclude) {
    # check to see if a combo contains an excluded pair
    # combo:
    #   - string name of combo (e.g. "BDQ+CLZ+JNJ")
    # pairs_exclude:
    #   - vector of pairs that should not be included (e.g. c("BDQ+JNJ", "TBA+TBI"))
    # returns TRUE if it does, FALSE if does not
    
    # break combo into pairs
    combo_singles <- combo_str2vec(combo)
    pairs_in_combo <- combos_to_text(combn(combo_singles, 2, simplify = FALSE))
    # check if any are in the pair_exclude list
    return(any(pairs_in_combo %in% pairs_exclude))
  } # end of combo_contains_excluded_pair
  
  combo_str2vec <- function(combo_text) {
    return(c(unlist(str_split(combo_text, "\\+"))))
  } # end of combo_str2vec
  
  combo_sort_drug_order <- function(combo) {
    # break apart combo and resort individual drug names
    return(paste(sort(combo_str2vec(combo)), collapse = "+"))
  } # end of combo_sort_drug_order
  
  get_combo_num_drugs <- function(combo) {
    return(length(combo_str2vec(combo)))
  } # end of get_combo_num_drugs
  
  
  # start main function
  # ============================================================================
  drug_list <- sort(drug_list) # sort into alphabetical order
  
  pairs_exclude <- c() # initialize empty vector to hold excluded pairs
  for (exclusion_category in drug_pairs_exclude) {
    excluded_pairs <- combn(sort(exclusion_category), 2, simplify = FALSE)
    pairs_exclude <- append(pairs_exclude, combos_to_text(excluded_pairs))
  }
  
  # how many n-ways of drugs to make
  combos_list <- list()
  for (n in nways) {
    combos_list <- append(combos_list, combn(drug_list, n, simplify = FALSE))
  }
  combos_list <- combos_to_text(combos_list)
  
  # scan through and exclude combos containing excluded pairs
  exclusion_list <- lapply(combos_list, combo_contains_excluded_pair, pairs_exclude)
  combos_list_filtered <- combos_list[exclusion_list == FALSE]
  
  # include labeled data
  labeled_combos <- labeled_data$Drug
  
  # double check the alphabetic sorting
  labeled_combos_sorted <- sapply(labeled_combos, combo_sort_drug_order, USE.NAMES = FALSE)
  
  combos_full <- rbind(
    tibble(combos_list_filtered) %>% rename(Drug = combos_list_filtered),
    tibble(labeled_combos_sorted) %>% rename(Drug = labeled_combos_sorted)) %>%
    unique()
  
  combos_full <- combos_full %>% rowwise() %>% mutate(NumbDrugs = get_combo_num_drugs(Drug))
  
  return(combos_full)
}


combos_to_text <- function(combos_list) {
  # applies a collapse function across a list of combos generated by combn to result in string
  # returns a vector with the names of the combos
  collapse_combo <- function(combo) {
    # takes c("BDQ", "CBR", "CLZ") and turns into "BDQ+CBR+CLZ"
    paste(combo, collapse = "+")
  }
  
  combos_list_text <- lapply(combos_list, collapse_combo)
  combos_list_text <- c(unlist(combos_list_text))
  return(combos_list_text)
} # end of combos_to_text


expand_combo2pairs <- function(combo_name) {
  # expands a string of high order combination name into a list of pairs
  # i.e. turns "A+B+C" into "A+B", "B+C", "A+C"
  # this is 5.8X faster than the old implementation
  drug_single_list <- c(unlist(str_split(combo_name, '\\+')))
  drug_pairs_comb <- combn(sort(drug_single_list), m = 2, simplify = FALSE)
  
  # function to generate pair list
  format_combn <- function(combn_element) {
    isolated_element = combn_element[[1]]
    concat_drugs <- paste(isolated_element, collapse = "+")
    return(concat_drugs)
  }
  
  pairs_list <- character(length(drug_pairs_comb)) # initialize
  
  # main loop
  for (i in 1:length(pairs_list)) {
    pairs_list[i] <- format_combn(drug_pairs_comb[i])
  }
  
  return(pairs_list)
}


aggregate_pairwise_diamond <- function(combo, diamond_pair_names, pairwise_metrics, metrics, n) {
  # make a vector that contains the data. MClapply will make it into lists
  # need to unpack the use of (i) as the function call
  # i in the main script is the Higher_order_comb (list of combo names)
  # inputs:
  #   - metrics
  #   - High_order_comb
  #   - diamond_pair_names
  #   - pairwise_metrics
  #   n: length of metrics

  # initialize an empty vector for performance
  combo_data <- rep(NA, n * 3)
  
  # expand out the combo to vector of pairs
  pairs_list <- expand_combo2pairs(combo)
  
  # check the match for each variable
  match_idx <- which(diamond_pair_names %in% pairs_list)
  matching_diamond_data <- pairwise_metrics[match_idx, ]
  
  # apply metrics calculations
  for (j in 1:n) {
    combo_data[j] <- max(matching_diamond_data[, j, drop = TRUE], na.rm = TRUE)
    combo_data[j+n] <- min(matching_diamond_data[, j, drop = TRUE], na.rm = TRUE)
    combo_data[j+n*2] <- mean(matching_diamond_data[, j, drop = TRUE], na.rm = TRUE)
  }
  return(combo_data)
}



aggregate_pairwise4combo <- function(combo_names, reformatted_diamond_data) {
    # This script takes the reformatted Diamond file and subsets higher order to lower order
    #Make sure to adjust for the input data points including df and f_name

    # combo_names actually has two columns: Drug and NumbDrugs
    stats_col <- c("Stdev","Average") # summary stats, don't need to split these for pairwise

    #extract pairwise drugs and important features
    pairwise <- reformatted_diamond_data %>% filter(NumbDrugs == 2) %>% select(-matches(paste(stats_col,collapse = "|")))
    # ^ pairwise data from DiaMond, same as above -stdev and average columns

    diamond_pair_names <- pairwise$Drug   # all pairwise drugs from compendium 
    pairwise_metrics <- pairwise %>% select(-c(Drug,NumbDrugs)) #everything except drugs and drug numb
    metrics <- names(pairwise_metrics) # characters of names of metrics

    # make column names with metrics
    col_names <- c(paste0(metrics, "_max"), 
                paste0(metrics, "_min"), 
                paste0(metrics, "_mean"))

    # parameters
    n <- length(metrics)

    # lapply the function
    pboptions(type="timer")
    b <- pblapply(c(combo_names$Drug), aggregate_pairwise_diamond, diamond_pair_names, pairwise_metrics, metrics, n, cl = detectCores()-2)
    b2 <- do.call(rbind, b) # combine all rows together
    b2[is.infinite(b2)] <- NA # delete data where unable to find max or min
    
    
    # assign column names and transform to a tibble for output
    colnames(b2) <- col_names
    b3 <- as_tibble(b2)
    
    
    # new function to count NumbDrugs in a text combo string
    num_drugs_from_name <- function(combo_name) {
        return(length(c(unlist(str_split(combo_name, "\\+")))))
    }
    
    num_drugs <- lapply(c(combo_names$Drug), num_drugs_from_name)
    b3 <- b3 %>% add_column(NumbDrugs = num_drugs, .before = 1)
    
    final_data <- b3 %>% add_column(Drug = c(combo_names$Drug), .before = 1)
    final_data$NumbDrugs <- c(unlist(final_data$NumbDrugs)) # to enable sorting
    return(final_data)
}


get_predictions <- function(combo_names, diamond_data_reformatted, subconditions, final_mod_list, Xproj_list, df_ThreshVPerf_list) {
    # should be able to get predictions from the name of the combo, reformatted diamond_data, and subcondition model

    # first, make a list of variables that are called from the environment
    # the pred_fun seems to contain everything
    # (cond_name, model_list, Xproj_list, ThreshVPerf_list, orig_pred_list, data_4_pred)
    # orig_pred_list seems not to be used

    # AUCandF1_list
    # bmr_class_lrn1

    #scripts: takes input the model training script workspace, and the list of combos that are labled and the list of 
    #combos that we need to make predition, and a list of models to use for making prediction 

    # rename function
    var_rename <- function(data) {
        data <- data %>%
            # replace condition
            rename_at(vars(matches("Butyrate")), list(~str_replace_all(., "Butyrate", "b"))) %>%
            rename_at(vars(matches("Valerate")), list(~str_replace_all(., "Valerate", "v"))) %>%
            rename_at(vars(matches("acidic")), list(~str_replace_all(., "acidic", "a"))) %>%
            rename_at(vars(contains("low")), list(~str_replace_all(., "lowcholesterol", "c"))) %>%
            rename_at(vars(matches("cholesterol")), list(~str_replace_all(., "cholesterol", "h"))) %>%
            rename_at(vars(matches("intracellular")), list(~str_replace_all(., "intracellular", "i"))) %>%
            rename_at(vars(matches("standard")), list(~str_replace_all(., "standard", "s"))) %>%
            rename_at(vars(matches("dormancy")), list(~str_replace_all(., "dormancy", "d"))) %>%
            # replace time points
            rename_at(vars(matches("Constant")), list(~str_replace_all(., "Constant", "C"))) %>%
            rename_at(vars(matches("ConstantTerminal")), list(~str_replace_all(., "ConstantTerminal", "CT"))) %>%
            rename_at(vars(matches("Terminal")), list(~str_replace_all(., "Terminal", "T")))
        return(data)
    }

    # prediction function
    pred_fun <- function(cond_name, model_list, Xproj_list, ThreshVPerf_list, data_4_pred) {
        # select subset data
        Xproj <- Xproj_list[[cond_name]] # original data
        mod<-model_list[[cond_name]] # model
        TvP<-ThreshVPerf_list[[cond_name]] # threshold vs performance 
        
        # training drugs and threshold performance data
        train_drug <- Xproj %>% filter(Dataset == "Train") %>% pull(Drug)
        ThreshVPerf <-df_ThreshVPerf_list[[cond_name]]$data # performance v threshold data
        do_pred <- predict( object = mod, newdata = data_4_pred %>% as.data.frame)
        pred_drugs <- data_4_pred %>% select(Drug, NumbDrugs) # the three way drugs you are testing on   ########
        
        # for Youden's J
        # # get averaged performance data from the resampling/crossvalidation and calculate youden's J
        train_perf_data <-TvP$data %>% # get resampling data
            mutate(Y_J = tpr+tnr-1) #  calculate youden's J for each threshold value
        
        # get optimal threshold from Youden's J
        all_Y_J_max <- train_perf_data%>%
            slice_max(Y_J) %>% pull(threshold) # select the threshold where Youden's J is maximum
        
        # optimum threshold
        opt_thresh <-all_Y_J_max[all_Y_J_max!=0] %>% min # in the case of ties in the Max Youden's J, exclude 0 threshold and select the minimum (least restrictive)
        
        #Add the Drug Column, dataset (train/test), score add Youden's J threshold, 
        drug_dataset_score <- pred_drugs %>% mutate(Dataset = case_when(Drug %in% train_drug~"Train",TRUE~"Test")) %>%
            left_join(Xproj %>% select(Drug,Score),by = "Drug")
        extract_pred_data <- do_pred$data %>%
            add_column(YJ = opt_thresh) 
        extract_pred_data<- extract_pred_data %>% setNames(paste0(cond_name, "_",names(.))) # append the condition subset name to the variable names
        extract_pred_data <- bind_cols(drug_dataset_score,extract_pred_data)
        #output
        return(extract_pred_data)
    }

    #add the candidate combination with DiaMOND data 
    # data_4_predict can be built by the function aggregate_pairwise4combo
    data_4_predict <- aggregate_pairwise4combo(combo_names, diamond_data_reformatted)

    data_4_predict <- var_rename(data_4_predict) %>% unique

    #select condition subset model
    top_model_predictions <- data_4_predict %>% select(Drug, NumbDrugs)

    for (subcondition in subconditions) {
        print(paste0("    - Subcondition ", subcondition))
        cond_sub <- paste0("CondSub_", subcondition)
        All_pred_truelab <- pred_fun(
            cond_sub,
            final_mod_list,
            Xproj_list,
            df_ThreshVPerf_list,
            data_4_predict)
        
        All_pred_truelab <- All_pred_truelab[order(All_pred_truelab$Drug), ]
        
        All_pred_truelab3 <- All_pred_truelab %>%
            select( Drug, paste0("CondSub_", subcondition, "_prob.1"))
        top_model_predictions <- merge(top_model_predictions, All_pred_truelab3, by = "Drug")
    }

    return(top_model_predictions)
}

# rename function
var_rename <- function(data) {
  data <- data %>%
    # replace condition
    rename_at(vars(matches("Butyrate")), list(~str_replace_all(., "Butyrate", "b"))) %>%
    rename_at(vars(matches("Valerate")), list(~str_replace_all(., "Valerate", "v"))) %>%
    rename_at(vars(matches("acidic")), list(~str_replace_all(., "acidic", "a"))) %>%
    rename_at(vars(contains("low")), list(~str_replace_all(., "lowcholesterol", "c"))) %>%
    rename_at(vars(matches("cholesterol")), list(~str_replace_all(., "cholesterol", "h"))) %>%
    rename_at(vars(matches("intracellular")), list(~str_replace_all(., "intracellular", "i"))) %>%
    rename_at(vars(matches("standard")), list(~str_replace_all(., "standard", "s"))) %>%
    rename_at(vars(matches("dormancy")), list(~str_replace_all(., "dormancy", "d"))) %>%
    # replace time points
    rename_at(vars(matches("Constant")), list(~str_replace_all(., "Constant", "C"))) %>%
    rename_at(vars(matches("ConstantTerminal")), list(~str_replace_all(., "ConstantTerminal", "CT"))) %>%
    rename_at(vars(matches("Terminal")), list(~str_replace_all(., "Terminal", "T")))
  return(data)
}





get_predictions_from_data <- function(high_order_data, subconditions, final_mod_list, Xproj_list, df_ThreshVPerf_list) {
  # get predictions from a higher-order dataset (data_4_predict)
  # first, make a list of variables that are called from the environment
  # the pred_fun seems to contain everything
  # (cond_name, model_list, Xproj_list, ThreshVPerf_list, orig_pred_list, data_4_pred)
  # orig_pred_list seems not to be used
  
  # AUCandF1_list
  # bmr_class_lrn1
  
  #scripts: takes input the model training script workspace, and the list of combos that are labled and the list of 
  #combos that we need to make predition, and a list of models to use for making prediction 
  
  
  
  # prediction function
  pred_fun <- function(cond_name, model_list, Xproj_list, ThreshVPerf_list, data_4_pred) {
    Rprof()
    # select subset data
    Xproj <- Xproj_list[[cond_name]] # original data
    mod <- model_list[[cond_name]] # model
    TvP <- ThreshVPerf_list[[cond_name]] # threshold vs performance 
    
    # training drugs and threshold performance data
    train_drug <- Xproj %>% filter(Dataset == "Train") %>% pull(Drug)
    ThreshVPerf <-df_ThreshVPerf_list[[cond_name]]$data # performance v threshold data
    do_pred <- predict(object = mod, newdata = data_4_pred %>% as.data.frame)
    pred_drugs <- data_4_pred %>% select(Drug, NumbDrugs) # the three way drugs you are testing on   ########
    
    # for Youden's J
    # # get averaged performance data from the resampling/crossvalidation and calculate youden's J
    train_perf_data <-TvP$data %>% # get resampling data
      mutate(Y_J = tpr+tnr-1) #  calculate youden's J for each threshold value
    
    # get optimal threshold from Youden's J
    all_Y_J_max <- train_perf_data %>%
      slice_max(Y_J) %>%
      pull(threshold) # select the threshold where Youden's J is maximum
    
    # optimum threshold
    opt_thresh <-all_Y_J_max[all_Y_J_max != 0] %>% min # in the case of ties in the Max Youden's J, exclude 0 threshold and select the minimum (least restrictive)
    
    #Add the Drug Column, dataset (train/test), score add Youden's J threshold, 
    drug_dataset_score <- pred_drugs %>%
      mutate(Dataset = case_when(Drug %in% train_drug ~ "Train",TRUE ~ "Test")) %>%
      left_join(Xproj %>% select(Drug,Score), by = "Drug")
    extract_pred_data <- do_pred$data %>%
      add_column(YJ = opt_thresh) 
    extract_pred_data<- extract_pred_data %>% setNames(paste0(cond_name, "_",names(.))) # append the condition subset name to the variable names
    extract_pred_data <- bind_cols(drug_dataset_score, extract_pred_data)
    #output
    Rprof(NULL)
    return(extract_pred_data)
  }
  
  #add the candidate combination with DiaMOND data 
  # data_4_predict can be built by the function aggregate_pairwise4combo
  data_4_predict <- var_rename(high_order_data) %>% unique
  
  #select condition subset model
  top_model_predictions <- data_4_predict %>% select(Drug, NumbDrugs)
  
  for (subcondition in subconditions) {
    print(paste0("    - Subcondition ", subcondition))
    cond_sub <- paste0("CondSub_", subcondition)
    All_pred_truelab <- pred_fun(
      cond_sub,
      final_mod_list,
      Xproj_list,
      df_ThreshVPerf_list,
      data_4_predict)
    
    All_pred_truelab <- All_pred_truelab[order(All_pred_truelab$Drug), ]
    
    All_pred_truelab3 <- All_pred_truelab %>%
      select( Drug, paste0("CondSub_", subcondition, "_prob.1"))
    top_model_predictions <- merge(top_model_predictions, All_pred_truelab3, by = "Drug")
  }
  
  return(top_model_predictions)
}

feat_imp_subcondition <- function(task_list, subcondition) {
  # needs the task_list from the model

  n = 100 # takes the average feature importance across 100 iterations
  condition_name_formatted <- paste0("CondSub_", tolower(subcondition))
  feat_imp <- replicate(n, generateFilterValuesData(task_list[[condition_name_formatted]], method = "randomForestSRC_importance"), simplify = FALSE)

  # grab data frames "data" from output list and rename "value" column
  data <- lapply(feat_imp, "[[", "data")
  data <- lapply(seq_along(data), function(x) setNames(data[[x]], c("name", "type", "filter", paste0("value", x))))

  # combine into one data frame and calculate median feature importance values
  meta_data <- data[[1]] %>% select(c("name", "type", "filter"))
  values <- bind_cols(lapply(data, function(x) x %>% select(starts_with("value"))))
  data <- bind_cols(list(meta_data, values))
  data$median_value <- apply(data %>% select(starts_with("value")), 1, median, na.rm = TRUE)
  data$sd_value <- apply(data %>% select(starts_with("value")), 1, sd, na.rm = TRUE)
  data$name <- reorder(data$name, desc(data$sd_value))

  return(data)
} # end of function


calculate_data_completeness <- function(combo_name, top_metrics, diamond_reformatted) {
  # calculates the data completeness score of a combo from pairwise diamond data
  # combo_name: name of a combo e.g. "BDQ+CLZ+PZA"
  # top_metrics: vector of metrics in shorthand form: Einf_b_C, GRinf_a_T
  # diamond_reformatted: diamond data reformatted to wide format
  expanded_combo <- expand_combo2pairs(test_combo)

  stats_col <- c("Stdev", "Average")

  pairwise <- diamond_reformatted %>%
    filter(NumbDrgs == 2) %>% 
    select(-matches(Paste(stats_col, collapse = "|")))

  # Have to group rename the indices
  pairwise <- rename_metric_to_abbreviated(pairwise)
}

rename_metric_to_abbreviated <- function(data) {
  # turns into shorthand e.g. Einf_Butyrate_Constant to Einf_b_C
  
  data <- data %>% 
    # replace condition
    rename_at(vars(matches("Butyrate")),list(~ str_replace(., "Butyrate", "b"))) %>% 
    rename_at(vars(matches("Valerate")),list(~ str_replace(., "Valerate", "v"))) %>% 
    rename_at(vars(matches("acidic")),list(~ str_replace(., "acidic", "a"))) %>% 
    rename_at(vars(contains("low")),list(~ str_replace(., "lowcholesterol", "c"))) %>% 
    rename_at(vars(matches("cholesterol")),list(~ str_replace(., "cholesterol", "h"))) %>% 
    rename_at(vars(matches("intracellular")),list(~ str_replace(., "intracellular", "i"))) %>% 
    rename_at(vars(matches("standard")),list(~ str_replace(., "standard", "s"))) %>% 
    rename_at(vars(matches("dormancy")),list(~ str_replace(., "dormancy", "d"))) %>%
    # replace time points
    rename_at(vars(matches("Constant")),list(~ str_replace(., "Constant", "C"))) %>%
    rename_at(vars(matches("ConstantTerminal")),list(~ str_replace(., "ConstantTerminal", "CT"))) %>%
    rename_at(vars(matches("Terminal")),list(~ str_replace(., "Terminal", "T")))
  return(data)
}

get_data_completeness <- function(combo_name, top_metrics, pairwise_diamond_data) {
  # select metrics and rows
  selected_metrics <- get_data_completeness_matrix(combo_name, top_metrics, pairwise_diamond_data)

  # data completeness score
  n_rows <- nrow(selected_metrics)
  n_cols <- ncol(selected_metrics) - 1 # because of the Drug column

  numel <- n_rows * n_cols # total number of elements
  n_na <- sum(is.na(selected_metrics))

  # have to account for elements that are from missing pairs
  # assume they are all na
  n_missing_pairs <- length(expand_combo2pairs(combo_name)) - n_rows
  n_na <- n_na + (n_missing_pairs * n_cols)
  numel <- numel + (n_missing_pairs * n_cols)

  data_completeness_score <- (numel - n_na) / numel
  return(data_completeness_score)
}

get_data_completeness_matrix <- function(combo_name, top_metrics, pairwise_diamond_data) {
  # calculate data completeness score for a specific combo
  expanded_combo <- expand_combo2pairs(combo_name)

  # select metrics and rows
  selected_metrics <- pairwise_diamond_data %>% select(c(Drug, all_of(top_metrics))) %>% filter(Drug %in% expanded_combo)
  return(selected_metrics)
}