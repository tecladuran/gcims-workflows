# baseline.R — GC-IMS Baseline Correction Functions

#' Add fixed-size area to the peak list
#'
#' Computes the fixed-size area of each peak based on its RT and DT range
#' and adds it as a new column `fixedsize_Area`.
#'
#' @param peak_list Data frame containing fixedsize_dt_* and fixedsize_rt_* columns
#'
#' @return Updated peak_list with fixedsize_Area column
#' @export
addFixedSizeArea <- function(peak_list) {
  peak_list <- peak_list %>%
    mutate(
      fixedsize_dt_range = fixedsize_dt_max_ms - fixedsize_dt_min_ms,
      fixedsize_rt_range = fixedsize_rt_max_s - fixedsize_rt_min_s,
      fixedsize_Area = fixedsize_dt_range * fixedsize_rt_range
    )

  return(peak_list)
}


#' Extract minimum intensity values for each peak region
#'
#' Computes the 5th percentile of intensity values within each peak region of the peak list.
#' The patch can be expanded by a percentage to better estimate the background.
#'
#' @param dataset GC-IMS dataset loaded with GCIMS package
#' @param peak_list Data frame with peak coordinates (must include SampleID and peak boundaries)
#' @param ampliation Percentage to expand the peak window (default = 0)
#'
#' @return Updated peak_list with a MinValue column
#' @export
extractMinValues <- function(dataset, peak_list, ampliation = 0) {
  peak_list$MinValue <- NA_real_

  for (sample_id in unique(peak_list$SampleID)) {
    current_sample <- dataset$getSample(sample_id)
    current_peaks <- dplyr::filter(peak_list, SampleID == sample_id)

    for (i in seq_len(nrow(current_peaks))) {
      dt_range <- c(current_peaks$fixedsize_dt_min_ms[i], current_peaks$fixedsize_dt_max_ms[i])
      rt_range <- c(current_peaks$fixedsize_rt_min_s[i], current_peaks$fixedsize_rt_max_s[i])

      dt_exp <- (dt_range[2] - dt_range[1]) * ampliation / 100
      rt_exp <- (rt_range[2] - rt_range[1]) * ampliation / 100

      dt_range <- c(dt_range[1] - dt_exp/2, dt_range[2] + dt_exp/2)
      rt_range <- c(rt_range[1] - rt_exp/2, rt_range[2] + rt_exp/2)

      patch <- intensity(current_sample, dt_range = dt_range, rt_range = rt_range)
      minValue <- quantile(patch, 0.05, na.rm = TRUE)
      peak_list$MinValue[which(peak_list$SampleID == sample_id)[i]] <- max(minValue, 0)
    }
  }

  return(peak_list)
}


#' Compute residual volume for all clusters and samples
#'
#' Calculates baseline area estimation per cluster and sample based on a minimum intensity
#' and region size. Uses cluster_stats as fallback if no peak is found.
#'
#' @param dataset GC-IMS dataset
#' @param peak_list Data frame of peaks with cluster and sample IDs
#' @param cluster_stats Data frame with fallback info per cluster
#' @param ampliation Expansion percentage (default = 0)
#'
#' @return Data frame (matrix-style) with residual volume estimates
#' @export
computeResidualVolume <- function(dataset, peak_list, cluster_stats, ampliation = 0) {
  cluster_ids <- unique(peak_list$cluster)
  sample_ids <- unique(peak_list$SampleID)
  residual_volume_table <- matrix(NA, nrow = length(cluster_ids), ncol = length(sample_ids),
                                  dimnames = list(cluster_ids, sample_ids))

  for (sample_id in sample_ids) {
    current_sample <- dataset$getSample(sample_id)
    current_peaks <- dplyr::filter(peak_list, SampleID == sample_id)

    for (current_cluster in cluster_ids) {
      cluster_peaks <- dplyr::filter(current_peaks, cluster == current_cluster)

      if (nrow(cluster_peaks) > 0) {
        dt_range <- c(cluster_peaks$fixedsize_dt_min_ms[1], cluster_peaks$fixedsize_dt_max_ms[1])
        rt_range <- c(cluster_peaks$fixedsize_rt_min_s[1], cluster_peaks$fixedsize_rt_max_s[1])
        area <- cluster_peaks$fixedsize_Area[1]
      } else {
        cluster_info <- dplyr::filter(cluster_stats, cluster == current_cluster)
        if (nrow(cluster_info) == 0) next
        dt_range <- c(cluster_info$dt_min_ms[1], cluster_info$dt_max_ms[1])
        rt_range <- c(cluster_info$rt_min_s[1], cluster_info$rt_max_s[1])
        area <- cluster_info$dt_length_ms[1] * cluster_info$rt_length_s[1]
      }

      dt_exp <- (dt_range[2] - dt_range[1]) * ampliation / 100
      rt_exp <- (rt_range[2] - rt_range[1]) * ampliation / 100
      dt_range <- c(dt_range[1] - dt_exp/2, dt_range[2] + dt_exp/2)
      rt_range <- c(rt_range[1] - rt_exp/2, rt_range[2] + rt_exp/2)

      patch <- intensity(current_sample, dt_range = dt_range, rt_range = rt_range)
      minValue <- ifelse(length(patch) == 0, NA, max(quantile(patch, 0.05, na.rm = TRUE), 0))
      residual_volume <- minValue * area
      residual_volume_table[as.character(current_cluster), as.character(sample_id)] <- residual_volume
    }
  }

  return(as.data.frame(residual_volume_table))
}


#' Apply full baseline correction pipeline to a GC-IMS dataset
#'
#' This function performs baseline correction by:
#' 1. Adding fixed-size area to the peak list
#' 2. Estimating baseline using intensity quantiles
#' 3. Computing residual volumes per peak region
#' 4. Imputing missing values in the peak table
#' 5. Subtracting baseline volume from the peak table
#'
#' @param dataset GCIMSDataset object
#' @param peak_list Peak list with fixedsize limits and cluster assignments
#' @param cluster_stats Data frame with cluster-level fallback info
#' @param ampliation Expansion percentage used to estimate the baseline (default = 200)
#'
#' @return Baseline-corrected peak table as a data frame
#' @export
correctBaseline <- function(dataset, peak_list, cluster_stats, ampliation = 200) {
  peak_list <- addFixedSizeArea(peak_list)
  peak_list <- extractMinValues(dataset, peak_list, ampliation = ampliation)

  residual_volume_df <- computeResidualVolume(
    dataset = dataset,
    peak_list = peak_list,
    cluster_stats = cluster_stats,
    ampliation = ampliation
  )

  imputed_matrix <- imputePeakTable(
    peakTable(peak_list, aggregate_conflicting_peaks = max)$peak_table_matrix,
    dataset,
    cluster_stats
  )

  corrected_table <- as.data.frame(imputed_matrix) - t(residual_volume_df)

  return(corrected_table)
}
