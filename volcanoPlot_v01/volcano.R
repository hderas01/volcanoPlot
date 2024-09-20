library(ggsignif)
library(ggrepel)
library(tidyverse)
#library(plotly)


is_combo_contain_drugpair <- function(combo_name, drug_pair) {
  combo_singles <- c(unlist(str_split(combo_name, "\\+")))
  drug1 <- drug_pair[[1]][1]
  drug2 <- drug_pair[[1]][2]
  
  # if combo singles contains both drugs, then assign value 1, else 0
  if (drug1 %in% combo_singles && drug2 %in% combo_singles) {
    return(1)
  } else {
    return(0)
  }
}

is_combo_contain_drug <- function(combo_name, drug) {
  combo_singles <- c(unlist(str_split(combo_name, "\\+")))
  idx <- which(drug %in% combo_singles)

  if (identical(idx, integer(0))) {
    return(0) # no match found
  } else {
    return(1) # match found
  }
}

add_col_contains_drug <- function(df, drug_name) {
  new_col <- lapply(df$Drug, is_combo_contain_drug, drug = drug_name)
  df <- add_column(df, new_col, .before = 2)
  df <- rename(df, !!paste("contains_", drug_name, sep = "") := new_col)
  return(df)
}

get_drug_pair_name <- function(drug_pair) {
  drug_pair_name <- paste0(drug_pair[[1]], "+", drug_pair[[2]])
  #print(drug_pair_name)
}


get_drug_pair_n <- function(drugpair) {
  drug_pair_col <- mclapply(df_drugpair$Drug, is_combo_contain_drugpair, 
                            drug_pair = drugpair,
                            mc.cores = detectCores()-2)
  num_drugpair_match <- sum(drug_pair_col == 1)
  return(num_drugpair_match)
}


# function that takes the main df and turns it into summary table
generate_drug_pair_table <- function(df_drugpair) {
  
  # inputs
  # df_drugpair
  # drugs
  
  # add column for drug_pair_col
  df_drugpair <- df_drugpair %>% add_column(drug_pair_col = 0)
  
  # make list of drug pairs
  drug_pairs_list <- combn(sort(unlist(drugs$abbreviation)), m = 2, simplify = FALSE)
  
  # filter out for only combos with 3 or 4 drugs
  df_drugpair <- df_drugpair %>% mutate(in_combo = if_else(drug_pair_col > 0, "Combos with pair", "Combos without pair"))
  
  # apply to all
  pair_name <- lapply(drug_pairs_list, get_drug_pair_name)
  drug_pair_table <- tibble(pair_name,
                            drug_pairs_list,
                            n = 0,
                            n_wo = 0,
                            prob_with = 0,
                            prob_wo = 0,
                            prob_diff = 0,
                            fold_change = 0,
                            P_wilcox = 0,
                            P_wilcox_log10 = 0)
  
  # main loop to calculate things
  for (i in 1:nrow(drug_pair_table)) {
    output <- get_drug_pair_metrics(drug_pair_table$drug_pairs_list[i], df_drugpair = df_drugpair)
    #print(output)
    drug_pair_table$n[i] = output[1]
    drug_pair_table$n_wo[i] = output[2]
    drug_pair_table$prob_with[i] = output[3]
    drug_pair_table$prob_wo[i] = output[4]
    drug_pair_table$prob_diff[i] = output[5]
    drug_pair_table$fold_change[i] = output[6]
    drug_pair_table$P_wilcox[i] = output[7]
    drug_pair_table$P_wilcox_log10[i] = output[8]
  }
  
  # remove rows that are 0, removes 10 entries
  drug_pair_table <- drug_pair_table %>% filter(n != 0)
  #drug_pair_table <- as_tibble(drug_pair_table)
  
  # apparently dataframe contains lists instead of vectors
  # conversion
  df2 <- drug_pair_table %>% select(-drug_pairs_list)
  df3 <- as_tibble(lapply(df2, unlist))
  
  # need a log2 fold change column
  drug_pair_table <- df3 %>% mutate(fold_change_log2 = log2(fold_change))
  return(drug_pair_table)
}

is_combo_contain_3_way <- function(combo_name, combo3) {

  combo_singles <- c(unlist(str_split(combo_name, "\\+")))
  drug1 <- combo3[[1]][1]
  drug2 <- combo3[[1]][2]
  drug3 <- combo3[[1]][3]
  
  # if combo singles contains 3 drugs, then assign value 1, else 0
  if (drug1 %in% combo_singles & drug2 %in% combo_singles & drug3 %in% combo_singles) {
    return(1)
  } else {
    return(0)
  }
}

is_combo_contain_subcombo <- function(combo, subcombo) {
  # break combination into singles
  combo_singles <- c(unlist(str_split(combo, "\\+")))
  subcombo_singles <- c(unlist(str_split(subcombo, "\\+")))

  if (all(subcombo_singles %in% combo_singles) == TRUE) {
    return(1)
  } else {
    return(0)
  }
}

is_combo_contain_any_subcombos <- function(combo, subcombo_list) {
  # break combination into singles
  combo_singles <- c(unlist(str_split(combo, "\\+")))

  # make subcombos from combo_singles
  subcombos_in_combo <- combos_to_text(combn(sort(combo_singles), m = 2, simplify = FALSE))

  return(sum(subcombos_in_combo %in% subcombo_list))
}

# main 3-way metrics script
get_combo3_metrics <- function(subcombo, df) {
  # setup list that shows location of where things are
  combo3_col <- mclapply(df$Drug, is_combo_contain_subcombo,
                            subcombo = subcombo,
                            mc.cores = detectCores()-2)
  df$combo3_col <- combo3_col
  

  # metrics
  n <- sum(combo3_col == 1)
  n_wo <- nrow(df) - n
  prob_with <- median(df$median_P[df$combo3_col == 1])
  prob_wo <- median(df$median_P[df$combo3_col == 0])
  prob_diff <- prob_with - prob_wo
  fold_change <- prob_with / prob_wo
  
  # statistics
  # only works if n > 0
  if (n > 0) {
    result <- wilcox.test(df$median_P[df$combo3_col == 1],
                        df$median_P[df$combo3_col == 0])
    P_wilcox <- result$p.value
    P_wilcox_log10 <- -log10(P_wilcox)
  } else {
    P_wilcox <- 0
    P_wilcox_log10 <- 0
  }
  
  
  output <- list(n,
                 n_wo, 
                 prob_with, 
                 prob_wo, 
                 prob_diff, 
                 fold_change, 
                 P_wilcox, 
                 P_wilcox_log10)
  return(output)
}

generate_n_way_table <- function(n, df) {
  # inputs
  # df_drugpair
  # drugs
  
  # add column for drug_pair_col
  #df <- df %>% add_column(combo3_present = 0)
  
  # make list of drug pairs THIS ALSO HAD TO BE CHANGED, DOUBLE CHECK
  singles_list <- get_singles_from_combos(df$Drug)
  subcombo <- combos_to_text(combn(sort(singles_list), m = n, simplify = FALSE))
  
  
  # recode
  #df <- df %>% mutate(in_combo = if_else(combo3_present > 0, "Combos with 3-way", "Combos without 3-way"))
  
  
  combo3_table <- tibble(subcombo,
                        n = 0,
                        n_wo = 0,
                        prob_with = 0,
                        prob_wo = 0,
                        prob_diff = 0,
                        fold_change = 0,
                        P_wilcox = 0,
                        P_wilcox_log10 = 0)
  
  # main loop to calculate things
  for (i in 1:nrow(combo3_table)) {
    output <- get_combo3_metrics(combo3_table$subcombo[i], df = df)
    combo3_table$n[i] = output[1]
    combo3_table$n_wo[i] = output[2]
    combo3_table$prob_with[i] = output[3]
    combo3_table$prob_wo[i] = output[4]
    combo3_table$prob_diff[i] = output[5]
    combo3_table$fold_change[i] = output[6]
    combo3_table$P_wilcox[i] = output[7]
    combo3_table$P_wilcox_log10[i] = output[8]
  }
  
  # remove rows that are 0, removes 10 entries
  combo3_table <- combo3_table %>% filter(n != 0)
  #drug_pair_table <- as_tibble(drug_pair_table)
  
  # apparently dataframe contains lists instead of vectors
  # conversion
  df2 <- combo3_table
  df3 <- as_tibble(lapply(df2, unlist))
  
  # need a log2 fold change column
  combo3_table <- df3 %>% mutate(fold_change_log2 = log2(fold_change))
  return(combo3_table)
  }



plot_pairs_volcano <- function(drug_pair_table, fc_thresh = 1.5, sig_thresh = 0.05) {
  
  volcano <- drug_pair_table
  
  # bin based on P-value first
  volcano <- volcano %>%
    mutate(signif_pair = cut(P_wilcox, breaks = c(0, sig_thresh, Inf)))
  
  # relabel to signif and not signif
  volcano <- volcano %>%
    mutate(signif_pair = case_match(signif_pair,
               paste0("(0,", sig_thresh, "]") ~ "Signif",
               paste0("(", sig_thresh, ",Inf]") ~ "Not signif"))
  
  volcano <- volcano %>%
    mutate(label = case_when(
      signif_pair == "Signif" & fold_change_log2 >= log2(fc_thresh) ~ "Outperforms",
      signif_pair == "Signif" & fold_change_log2 < log2(1/fc_thresh) ~ "Underperforms",
      .default = "Not significant"
    ))
  
  # label if it meets the threshold
  volcano$pair_label <- NA
  volcano$pair_label[volcano$label != "Not significant"] <- 
    volcano$pair_name[volcano$label != "Not significant"]
  
  # alternative: label everything
  volcano$pair_label <- volcano$pair_name
    
  # Plot volcano plot
  mycolors <- c("blue", "red", "darkgrey")
  names(mycolors) <- c("Underperforms", "Outperforms", "Not significant")
  
  p <- ggplot(volcano, 
              aes(
                x = fold_change_log2, 
                y = P_wilcox_log10, 
                color = label, 
                label = pair_label)) +
    geom_point() +
    theme_classic() +
    #geom_vline(xintercept = c(log2(1/fc_thresh), log2(fc_thresh)), col = "brown4", linewidth = 0.2) +
    geom_hline(yintercept = -log10(sig_thresh), col = "brown4", linewidth = 0.2) + 
    scale_color_manual(values = mycolors) +
    geom_text_repel(
      size = 2.5,
      max.overlaps = 18
    ) + 
    #xlim(-1.5, 1.0) + 
    #ylim(-10, 120) + 
    xlab("Log2(Fold change)") + 
    ylab("Log10(P value)")
  print(p)
  return(p)
}


plot_volcano_GC <- function(drug_pair_table, fc_thresh = 1.5, sig_thresh = 0.05) {
  
  volcano <- drug_pair_table
  
  # bin based on P-value first
  volcano <- volcano %>%
    mutate(signif_pair = cut(P_wilcox, breaks = c(0, sig_thresh, Inf)))
  
  # relabel to signif and not signif
  volcano <- volcano %>%
    mutate(signif_pair = case_match(signif_pair,
               paste0("(0,", sig_thresh, "]") ~ "Signif",
               paste0("(", sig_thresh, ",Inf]") ~ "Not signif"))
  
  volcano <- volcano %>%
    mutate(label = case_when(
      signif_pair == "Signif" & fold_change_log2 >= log2(fc_thresh) ~ "Outperforms",
      signif_pair == "Signif" & fold_change_log2 < log2(1/fc_thresh) ~ "Underperforms",
      .default = "Not significant"
    ))
  
  # label if it meets the threshold
  volcano$pair_label <- NA
  volcano$pair_label[volcano$label != "Not significant"] <- 
  volcano$pair_name[volcano$label != "Not significant"]
  
  # alternative: label everything
  #volcano <- volcano |> mutate(pair_label = case_when(
  #  pair_name == "JNJ+LIN+MMV" ~ pair_name,
  #  pair_name == "DEL+MMV+MOX" ~ pair_name,
  #  pair_name == "DEL+MMV+SFN" ~ pair_name,
  #  pair_name == "LIN+MMV+SFN" ~ pair_name,
  #  .default = NA
  #))
    
  # Plot volcano plot
  mycolors <- c("grey", "blue")
  #names(mycolors) <- c("Underperforms", "Outperforms", "Not significant")
  
  p <- ggplot(volcano, 
              aes(
                x = fold_change_log2, 
                y = P_wilcox_log10, 
                color = drug_int, 
                #changed color line from "contains_MOX" to "drug_int" 18Dec2023
                label = pair_label)) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_vline(xintercept = c(log2(1/fc_thresh), log2(fc_thresh)), col = "brown4", linewidth = 0.2) +
    geom_hline(yintercept = -log10(sig_thresh), col = "brown4", linewidth = 0.2) + 
    scale_color_manual(values = mycolors) +
    #geom_text_repel(
     # size = 3,
      #max.overlaps = 18,
      #color = "black",
      #arrow = arrow(length = unit(0.02, "npc")),
      #nudge_x = 0.5
    #) + 
    #xlim(-1.5, 1.0) + 
    #ylim(-10, 120) + 
    xlab("Log2(Fold change)") + 
    ylab("Log10(P value)")
  print(p)
  plotly(p)
  return(p)
}



# create list of all possible drug pairs