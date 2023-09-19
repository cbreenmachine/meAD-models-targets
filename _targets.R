# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    # argparse
    "DSS",
    "magrittr",
    "dplyr",
    "stringr",
    "parallel",
    "tidyverse",
    "cowplot",
    "qs",
    "data.table",
    "magrittr",
    "crew"
  ), 
  format = "qs",
  # https://books.ropensci.org/targets/crew.html
  controller = crew::crew_controller_local(workers = 4),
)

# tar_make_clustermq() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
options(clustermq.scheduler = "multicore")

# tar_make_future() is an older (pre-{crew}) way to do distributed computing
# in {targets}, and its configuration for your machine is below.
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

tar_source("R/")

# Directories and file paths
data_dir <- "data"
m_cov_dir <- file.path(data_dir, "methy-coverage/")
samplesheet <- file.path(data_dir, "masterSamplesheet.csv")

# Constants
chrom_df <- data.frame(chrom = paste0("chr", 1:22))

# Replace the target list below with your own:
prelims <- list(
  tar_target(master_df, clean_master_samplesheet("DataRaw/masterSamplesheet.csv")),
  tar_target(ids, get_valid_ids(master_df, "DataRaw/MCovRawSplit/")),
  tar_target(control_ids, get_control_samples(master_df, ids)),
  tar_target(load_ids, get_load_samples(master_df, ids))
)


data_filtering <- tar_map(
  values = chroms_df,
  tar_target(m_cov, 
             read_and_pivot_routine(m_cov_dir, chrom, ids)),
  tar_target(loci_stats,
             compute_loci_summary_stats(m_cov)),
  tar_target(by_sample_coverage_stats,
              compute_coverage_stats(m_cov)),
  tar_target(loci_info,
             generate_filter_ixs(m_cov, loci_stats)),
  tar_target(m_cov_imputed_filtered,
             impute_routine(m_cov, control_ids, load_ids, loci_info))
)

list(
  prelims,
  data_filtering
)
