# _targets.R
# Air traffic controller of pipeline used for meAD-paper
library(targets)
library(tarchetypes)

########################################################
###################### Parameters ######################
########################################################
LFDR.CUT <- 0.05

# Covariates to adjust for
baseline_covariates <- c("age", "is_male", "BMI", "race_primary")
derived_covariates <- c("PC1", "PC2", "type1", "type2", "type3", "type4", "type5")

# Testing (LOAD related)
TEST_COVARIATE <- "diagnostic_group_coded"

controller_high_mem <- crew::crew_controller_local(
    name = "high_mem_controller",
    workers = 1, 
    # BCG file system is slow, so if errors thrown increase this
    seconds_idle = 30
)

controller_default <- crew::crew_controller_local(
  name = "default_controller"
  workers = 4,
  seconds_idle = 10
)

########################################################
######################## Targets #######################
########################################################
tar_option_set(
  packages = c(
    "annotatr",
    "chunkR",
    "cowplot",
    "crew",
    "data.table",
    "devEMF",
    "dplyr",
    "DSS",
    "EnsDb.Hsapiens.v86",
    "fastqq",
    "GenomicRanges",
    "ggsci",
    "gridExtra",
    "harmonicmeanp",
    "htmlwidgets",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "minfi",
    "magrittr",
    "networkD3",
    "openxlsx",
    "parallel",
    "qs",
    "readxl",
    "rjson",
    "stringr",
    "tidyverse",    
    "tibble",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "viridis",
    "waffle"
  ), 
  format = "qs",
  controller = crew_controller_group(
    controller_default,
    controller_high_mem
  ),
  resources = tar_resources(
    crew = tar_resources_crew(controller = "default_controller")
  )
)



tar_source("R/")

# Directories and file paths
data_dir <- "data"
ref_dir <- "data_reference"
m_cov_dir <- file.path(data_dir, "methy_coverage/")
samplesheet <- file.path(data_dir, "master_samplesheet.csv")

# Construct a data frame to loop over in parallel
# so we can run multiple targets at once
#TODO: change this to sweep all chromosomes
chrom_df <- data.frame(chrom = paste0("chr", 1:22))

# Replace the target list below with your own:
prelims <- list(
  tar_target(master_full_df, prepare_samplesheet(samplesheet)),
  tar_target(master_df, clean_master_samplesheet(master_full_df, baseline_covariates)),
  tar_target(ids, get_valid_ids(master_df, m_cov_dir)),
  tar_target(control_ids, get_control_samples(master_df, ids)),
  tar_target(load_ids, get_load_samples(master_df, ids))
)


test_demographics <- list(
  
  tar_target(ancestry_test, tabulate_and_test(master_full_df, "race_primary", "fisher")),
  tar_target(source_test, tabulate_and_test(master_full_df, "source", "chi")),
  tar_target(sex_test, tabulate_and_test(master_full_df, "sex", "chi")),
  tar_target(APOE_test, tabulate_and_test(master_full_df, "APOE_risk_allele", "chi")),

  tar_target(age_test, test_continuous(master_full_df, "age_at_visit")),
  tar_target(bmi_test, test_continuous(master_full_df, "bmi")),
  tar_target(education_test, test_continuous(master_full_df, "education"))
)


data_filtering <- tar_map(
  values = chrom_df,

  # A length-two list with first item representing the methylated counts
  # and second item representing the total counts
  tar_target(m_cov, 
             read_and_pivot_routine(m_cov_dir, chrom, ids),
             resources = tar_resources(
              crew = tar_resources_crew(controller = "high_mem_controller")
    )
  ),
  
  # A length-two list with percent nonzero and median coverage
  tar_target(loci_stats,
             compute_loci_summary_stats(m_cov)),

  # Curiosity--what is the min, med, max coverage by sample
  tar_target(by_sample_coverage_stats,
              compute_coverage_stats(m_cov)),
  
  # Determines which loci are/are not filtered
  # First entry in list is the "keep index" and
  # second entry gives the percent of samples that pass each filter
  tar_target(loci_info,
             generate_filter_ixs(m_cov, loci_stats, percent_nonzero_thresh=0.5, med_cov_thresh=5)),

  # Given M/Cov matrices, and filter indices,
  # drop "bad" loci and impute as needed.
  # Imputa
  tar_target(m_cov_imputed_filtered,
             impute_routine(m_cov, control_ids, load_ids, loci_info)),
  
  tar_target(pcs, compute_pcs(m_cov_imputed_filtered))
)


# modeling <- list(
  # Make chromosome-specific design matrix (each gets a chrom-specific PC)

  # 
# )


# Reporting DMPs, effect sizes...
pvals_summary <- list(

)



list(
  prelims,
  test_demographics,
  data_filtering 
)
