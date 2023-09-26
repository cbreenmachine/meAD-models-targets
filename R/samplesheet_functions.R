


prepare_samplesheet <- function(ifile){
    # Remove extranesous columns, encode binary variables (sex, load status)
    # and check that there are no duplicate sample rows
    out_df <- read_csv(ifile, show_col_types = FALSE) %>%
        group_by(sample_id) %>%
        dplyr::slice(1) %>% 
        ungroup() %>%
        dplyr::arrange(sample_id) %>% 
        dplyr::mutate(sample_id = as.character(sample_id)) %>%
        dplyr::filter(diagnostic_group != "MCI") %>%
        dplyr::mutate(
            race_primary = as.factor(race_primary),
            APOE_risk_allele = (apoe_e1 == 4) + ( apoe_e2 == 4)
        ) %>%
        column_to_rownames("sample_id")

    if (length(unique(rownames(out_df))) == nrow(out_df)) {
        return(out_df)
    } else {
        stop("Incorrect number of samples in master_df")
    }
}


clean_master_samplesheet <- function(master_full_df, covariates){
    # Remove extranesous columns, encode binary variables (sex, load status)
    # and check that there are no duplicate sample rows
    out_df <- master_full_df %>%
        dplyr::mutate(
            age = age_at_visit,
            BMI = bmi,
            # Make CONTROL the reference level (0) to give us "LOAD effect"
            diagnostic_group_coded = ifelse(diagnostic_group == "CONTROL", 0, 1),
            # Female is the reference level (0) so this now gives "male effect"
            is_male = ifelse(sex == "Male", 1, 0),
            # Just a factor variable here
            race_primary = as.factor(race_primary)
        ) %>% 
        dplyr::select(all_of(c(covariates, "diagnostic_group_coded")))
        
    return(out_df)
}




tabulate_and_test <- function(master_df, var, test.type = "chi"){
  tb <- master_df %>%
    dplyr::count(diagnostic_group, !!rlang::sym(var)) %>%
    pivot_wider(names_from = "diagnostic_group", values_from = "n", values_fill = 0) %>%
    drop_na() %>%
    column_to_rownames(var)

  if (test.type == "chi"){
    test.result <- chisq.test(tb)
  } else if (test.type == "fisher"){
    test.result <- fisher.test(tb)
  }

  return(list("Table" = tb,
              "Test" = test.result))
}


get_valid_ids <- function(master_df, bed_dir){
    # Check that the samples in the directory have phenotype and vice versa
    idir_samples <- unique(stringr::str_split_fixed(list.files(bed_dir), "\\.", 3)[ ,2])

    # All samples that are in the directory and sampleshseet
    valid_samples <- intersect(rownames(master_df), idir_samples)
    valid_samples
}


get_load_samples <- function(master_df, valid_samples){
    load_samples <- intersect(
        rownames(master_df)[master_df$diagnostic_group == 1], 
        valid_samples)
    return(load_samples)
}

get_control_samples <- function(master_df, valid_samples){
    ctrl_samples <- intersect(
        rownames(master_df)[master_df$diagnostic_group_coded == 0], 
        valid_samples
    )
    return(ctrl_samples)
}



stringify_mean_sd <- function(xx){
    # Cleans up test outputs, could use broom:: in future versions
    paste0(round(mean(xx), 2), " (", round(sd(xx), 2), ")")
}


test_continuous <- function(samples_df, var){

  if (!(var %in% colnames(samples_df))){
    warning("Var not in df")
  }

  xx.control <- samples_df %>%
    dplyr::filter(diagnostic_group == "CONTROL") %>%
    dplyr::pull(var) %>% as.numeric()

  xx.load <- samples_df %>%
    dplyr::filter(diagnostic_group == "LOAD") %>%
    dplyr::pull(var) %>% as.numeric()

  list("Control: " = stringify_mean_sd(xx.control),
       "LOAD: " = stringify_mean_sd(xx.load),
       "Test: " = t.test(xx.control, xx.load))
}


