# # model_functions.R
# # Handles converting M/Cov matrices to BSSeq object
# # Fits DSS models

# ###########################################################
# ##################### DATA PREP ###########################
# ###########################################################
# covariates_split <- stringr::str_split(args$baseline_covariates, ",")[[1]]
# dss.formula <- formula(paste0("~", args$test_covariate, "+", paste(covariates_split, collapse = "+")))



# # Pull out just the covariates we'll use
# design.df <- df %>% 
#     tibble::column_to_rownames("sample_id") %>%
#     dplyr::select(all_of(c(args$test_covariate, covariates_split))) %>%
#     dplyr::mutate_at(args$test_covariate, as.numeric) %>% 
#     tidyr::drop_na()

# # If there are NAs in any of the covariates specified, be sure to drop
# #corresponding sample in bs object
# ix <- colnames(bs) %in% rownames(design.df)
# bs.sub <- bs[ , ix]

# if (all(colnames(bs.sub) == rownames(design.df))){
#     print("After filtering NAs, rownames in design match colnames in bs")
# } else {
#     warning("rownames in design != colnames in bs")
# }



# # Fit models
# dml.fit <- DMLfit.multiFactor(bs.sub, design = design.df, smoothing = smooth, 
#     smoothing.span = args$smoothing, formula = dss.formula)

# invisible(gc())

# # Test covariate (to start this will be LOAD status--diagnostic_group_coded)
# test.var <- args$test_covariate
# test.result <- DMLtest.multiFactor(dml.fit, coef = test.var)

# #To measure effect size, other beta coefficients...
# get_beta_estimate <- function(dml_fit){
#     return(data.frame(dml_fit$fit$beta))
# }

# beta.df <- 


# design.mat <- model.matrix(dss.formula, design.df)


# #END