#' Extract minimum intensity values for samples in a single cluster (exploratory use)
#'
#' Returns a vector of MinValues (5th percentile) for each sample in the specified cluster.
#' Useful for visualizing how baseline intensity changes with ampliation.
#'
#' @param dataset GC-IMS dataset
#' @param peak_list Data frame of peaks
#' @param cluster_name Cluster to analyze
#' @param ampliation Expansion percentage (default = 0)
#'
#' @return Numeric vector with minimum values per sample
#' @export
extractMinValues_2 <- function(dataset, peak_list, cluster_name, ampliation = 0) {
  cluster_peaks <- peak_list %>% filter(cluster == cluster_name)
  sample_ids <- unique(cluster_peaks$SampleID)
  min_values <- numeric(length(sample_ids))

  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[i]
    current_sample <- dataset$getSample(sample_id)

    row <- cluster_peaks[cluster_peaks$SampleID == sample_id, ][1, ]
    dt_range <- c(row$fixedsize_dt_min_ms, row$fixedsize_dt_max_ms)
    rt_range <- c(row$fixedsize_rt_min_s, row$fixedsize_rt_max_s)

    dt_exp <- (dt_range[2] - dt_range[1]) * ampliation / 100
    rt_exp <- (rt_range[2] - rt_range[1]) * ampliation / 100

    dt_range <- c(dt_range[1] - dt_exp/2, dt_range[2] + dt_exp/2)
    rt_range <- c(rt_range[1] - rt_exp/2, rt_range[2] + rt_exp/2)

    patch <- intensity(current_sample, dt_range = dt_range, rt_range = rt_range)
    min_values[i] <- ifelse(length(patch) == 0, NA, max(quantile(patch, 0.05, na.rm = TRUE), 0))
  }

  return(min_values)
}
