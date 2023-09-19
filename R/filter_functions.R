

.read_bed <- function(idir, chr, sample){
    #Creates input name like Data/chr22.100.bed
    # and reads it
    dt <- fread(file.path(idir, paste0(chr, ".", sample ,".bed")))
    dt$sample <- as.numeric(sample) # numeric quiets warning
    dt
}

pivot_me <- function(data, value) {
    # Gets into Sample \times Position matrix of 
    # M or Cov
    # "value" is "methylated" or "coverage"
    keepcols <- c("chromStart", "sample", value)

    data %>% 
        dplyr::select(all_of(keepcols)) %>%
        tidyr::pivot_wider(values_from = value, names_from = sample) %>%
        tibble::column_to_rownames("chromStart")
}


read_and_pivot_routine <- function(idir, chr, samples){
    # Read and stack all sample-specific data files
    .read_bed_wrapper <- function(sample){
        .read_bed(idir, chr, sample)
    }
    data <- do.call(rbind, lapply(X = samples, FUN = .read_bed_wrapper))

    m <- pivot_me(data, "methylated")
    cov <- pivot_me(data, "coverage")

    return(list(m = m, cov = cov))
}

compute_loci_summary_stats <- function(m_cov){
    m <- m_cov[["m"]]
    cov <- m_cov[["cov"]]

    # Copying neccessary inside the function scope?
    cov_zeroed <- data.table::copy(cov) 
    cov_zeroed[is.na(cov_zeroed)] <- 0

    # We'll use the number of samples as a divisor
    n_samples <- ncol(cov_zeroed)
    
    percent_zero <- rowSums(cov_zeroed == 0) / n_samples
    percent_nonzero <- 1 - percent_zero

    # median coverage calculation happens independently 
    # from zero
    med_cov <- apply(cov_zeroed, FUN = median, MARGIN = 1)

    return(list(percent_nonzero = percent_nonzero, 
                med_coverage = med_cov))
}


generate_filter_ixs <- function(m_cov, loci_stats, percent_nonzero_thresh=0.5, med_cov_thresh=5){
    m <- m_cov[["m"]]
    cov <- m_cov[["cov"]]

    # More unpacking
    percent_nonzero <- loci_stats[["percent_nonzero"]]
    med_coverage <- loci_stats[["med_coverage"]]

    # Number of CpGs
    G <- length(percent_nonzero)

    # non-zero and median coverage indices
    nz_ix <- (percent_nonzero >= percent_nonzero_thresh)
    mc_ix <- (med_coverage >= med_cov_thresh)

    keepix <- nz_ix & mc_ix

    stats <- data.frame(stat = c("perc_nonzero", "perc_ge_med_cov"),
                        value = c(nz_ix / G, mc_ix / G))
    
    return(list(keepix = keepix, stats = stats))
}


.filter_m_cov <- function(m_cov, keepix){
    m <- m_cov[["m"]]
    cov <- m_cov[["cov"]]

    m_filt <- data.table::copy(m)[keepix, ]
    cov_filt <- data.table::copy(cov)[keepix, ]

    return(list(m = m_filt, cov = cov_filt))
}


.create_filler_mask <- function(DT, ctrl.samples, load.samples, N){
    # Given the groups of samples, compute the mean for columns in 
    # "load.samples" and then make a data table that looks like DT,
    # but with columns imputed with mean. DT can be M or Cov

    # N is number of columns (samples)
    row.means.load <- round(rowMeans(DT[ , load.samples], na.rm = T))
    row.means.ctrl <- round(rowMeans(DT[ , ctrl.samples], na.rm = T))

    # Expand to make it a matrix we can mask
    filler <- DT
    
    filler[ , load.samples] <- row.means.load
    filler[ , ctrl.samples] <- row.means.ctrl

    return(filler)
}


impute_routine <- function(m_cov, ctrl.samples, load.samples, loci_info){

    keepix <- loci_info[["keepix"]]

    m_cov_filt <- .filter_m_cov(m_cov, keepix)

    m_filt <- m_cov_filt[['m']]
    cov_filt <- m_cov_filt[['cov']]

    # LOAD/Control row means expanded to be the same shape as M 
    m_means <- .create_filler_mask(m_filt, ctrl.samples, load.samples, N)
    cov_means <- .create_filler_mask(cov_filt, ctrl.samples, load.samples, N)

    # Masks for indexing
    # Should be the case that the values we want to impute are all zero by now
    # but include NA just in case
    mask <- (cov_filt == 0 | is.na(cov_filt))

    # Fill withe means
    m_filt[mask] <- m_means[mask]
    cov_filt[mask] <- cov_means[mask]

    return(list(m = m_filt, cov = cov_filt))
}


compute_coverage_stats <- function(m_cov){
    X <- m_cov[['cov']]

    mins <- apply(X, FUN=min, MARGIN=1)
    maxes <- apply(X, FUN=max, MARGIN=1)
    meds <- apply(X, FUN=median, MARGIN=1)
}
#END