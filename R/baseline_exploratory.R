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

extractMinValues_2 <- function(dataset, peak_list, cluster_name, ampliation = 0, n_samples) {

  # Filtrar només les entrades del cluster seleccionat
  cluster_peaks <- peak_list %>% filter(cluster == cluster_name)

  # Vector per emmagatzemar els valors mínims
  min_values <- numeric(n_samples)

  # Iterar sobre les mostres del cluster
  for (i in seq_len(n_samples)) {
    sample_id <- cluster_peaks$SampleID[i]
    # Obtenir la mostra corresponent
    current_sample <- dataset$getSample(sample_id)

    # Calcular els rangs de temps de deriva i de retenció
    dt_range <- c(cluster_peaks$fixedsize_dt_min_ms[i], cluster_peaks$fixedsize_dt_max_ms[i])
    rt_range <- c(cluster_peaks$fixedsize_rt_min_s[i], cluster_peaks$fixedsize_rt_max_s[i])

    # Calcular l'expansió
    dt_expansion <- (dt_range[2] - dt_range[1]) * ampliation / 100
    rt_expansion <- (rt_range[2] - rt_range[1]) * ampliation / 100

    # Aplicar l'expansió
    dt_range <- c(dt_range[1] - dt_expansion/2, dt_range[2] + dt_expansion/2)
    rt_range <- c(rt_range[1] - rt_expansion/2, rt_range[2] + rt_expansion/2)

    # Obtenir el fragment d'intensitat i calcular la mitjana del decil inferior
    patch <- intensity(current_sample, dt_range = dt_range, rt_range = rt_range)
    if (length(patch) == 0) {
      min_values[i] <- NA
    } else {
      min_values[i] <- quantile(patch, 0.05, na.rm = TRUE)
      min_values[i] <- max(min_values[i], 0)  # Ensure minValue is non-negative
    }
  }

  return(min_values)
}
