library(ggplot2)
library(cowplot)
library(gridExtra)

###########################################
#
# Build dataframe with clusters information
#
# parameters:
#  - peak_table: dataframe (rows: clusters, cols: samples)
# returns:
#  - clusters_info_df: dataframe with different characteristics for each cluster
###########################################
clusters_info <- function(peak_table){
  # Calculate number smples per cluster
  peak_table$samples <- rowSums(!is.na(peak_table[,-1]))
  clusters_info_df <- peak_table %>% select(cluster, samples)
  peak_table <- peak_table %>% select(-samples)

  # Calculate cv
  s <- apply(peak_table[,-1], 1, sd, na.rm=TRUE)
  m <- apply(peak_table[,-1], 1, mean, na.rm=TRUE)
  cv <- s/m
  clusters_info_df$cv <- cv

  # Calculate avg intenisty per cluster
  clusters_info_df$avg_int <- rowMeans(peak_table[,-1], na.rm = TRUE)
  return(clusters_info_df)
}

###########################################
#
# Minmax intensity to plot samples in same scale
#
# parameters:
#  - dataset
#  - idx_samples
# returns:
#  - (min_val, max_val)
###########################################
return_minmax_samples <- function(dataset, idx_samples){
  minims <- c()
  maxims <- c()

  trans <- cubic_root_trans()

  for (i in idx_samples){
    sample1 <- dataset$getSample(sample = i)
    intmat <- sample1@data
    intmat_trans <- trans$transform(intmat)
    minims <- c(minims, min(intmat_trans))
    maxims <- c(maxims, max(intmat_trans))
  }

  min_val <- min(minims)
  max_val <- max(maxims)
  return(c(min_val, max_val))
}

return_minmax_notransf_samples <- function(dataset, idx_samples){
  minims <- c()
  maxims <- c()

  trans <- cubic_root_trans()

  for (i in idx_samples){
    sample1 <- dataset$getSample(sample = i)
    intmat <- sample1@data
    minims <- c(minims, min(intmat))
    maxims <- c(maxims, max(intmat))
  }

  min_val <- min(minims)
  max_val <- max(maxims)
  return(c(min_val, max_val))
}

build_nr <- function(x, minmax = NULL){
  colormap <- farver::encode_native(
    viridisLite::viridis(256L, direction = -1, option = "G")
  )

  if (is.null(minmax)) {
    minmax <- range(x)
  }
  breaks <- seq(from = minmax[1], to = minmax[2], length.out = length(colormap))
  xdim <- dim(x)
  rev_cols <- seq.int(ncol(x), 1L, by = -1L)
  x <- x[, rev_cols]
  x <- findInterval(x, breaks, rightmost.closed = TRUE)
  x <- colormap[x]
  x <- matrix(x, nrow = xdim[1], ncol = xdim[2], byrow = FALSE)

  nr <- structure(
    x,
    dim = c(xdim[2], xdim[1]),
    class = "nativeRaster",
    channels = 4L
  )
  return(nr)
}

plot_sample <- function(nr, intmat, minmax, dt_min = NULL, dt_max = NULL, rt_min = NULL, rt_max = NULL, trans = "cubic_root"){
  if (is.null(dt_min)) {
    dt_min <- as.numeric(rownames(intmat)[1L])
  }
  if (is.null(dt_max)) {
    dt_max <- as.numeric(rownames(intmat)[nrow(intmat)])
  }
  if (is.null(rt_min)) {
    rt_min <- as.numeric(colnames(intmat)[1L])
  }
  if (is.null(rt_max)) {
    rt_max <- as.numeric(colnames(intmat)[ncol(intmat)])
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      xmin = dt_min, xmax = dt_min,
      ymin = rt_min, ymax = rt_min,
      ggplot2::aes(fill = .data$x),
      data = data.frame(
        x = NA_real_,
        dt_ms_min = dt_min, dt_ms_max = dt_max,
        rt_s_min = rt_min, rt_s_max = rt_max
      )
    ) +
    ggplot2::annotation_raster(
      nr,
      xmin = dt_min, xmax = dt_max,
      ymin = rt_min, ymax = rt_max
    ) +
    ggplot2::scale_fill_viridis_c( # This has to match with the COLORMAP above
      direction = -1,
      option = "G",
      limits = minmax,
      na.value = "#00000000",
      trans = trans
    ) +
    ggplot2::lims(
      x = c(dt_min, dt_max),
      y = c(rt_min, rt_max)
    ) +
    ggplot2::labs(
      x = "Drift time (ms)",
      y = "Retention time (s)",
      fill = "Intensity (a.u.)"
    ) +
    ggplot2::theme_minimal()
  print(p)
}

build_nr <- function(x, minmax = NULL, colormap){
  if (is.null(minmax)) {
    minmax <- range(x)
  }
  breaks <- seq(from = minmax[1], to = minmax[2], length.out = length(colormap))
  xdim <- dim(x)
  rev_cols <- seq.int(ncol(x), 1L, by = -1L)
  x <- x[, rev_cols]
  x <- findInterval(x, breaks, rightmost.closed = TRUE)
  x <- colormap[x]
  x <- matrix(x, nrow = xdim[1], ncol = xdim[2], byrow = FALSE)

  nr <- structure(
    x,
    dim = c(xdim[2], xdim[1]),
    class = "nativeRaster",
    channels = 4L
  )
  return(nr)
}

plot_sample <- function(sample_obj, colormap, minmax = NULL, minmax_for_legend = NULL, dt_min = NULL, dt_max = NULL, rt_min = NULL, rt_max = NULL){
  intmat <- intensity(sample_obj)

  if (is.null(dt_min)) {
    dt_min <- as.numeric(rownames(intmat)[1L])
  }
  if (is.null(dt_max)) {
    dt_max <- as.numeric(rownames(intmat)[nrow(intmat)])
  }
  if (is.null(rt_min)) {
    rt_min <- as.numeric(colnames(intmat)[1L])
  }
  if (is.null(rt_max)) {
    rt_max <- as.numeric(colnames(intmat)[ncol(intmat)])
  }

  trans <- cubic_root_trans()
  intmat_trans <- trans$transform(intmat)
  if (is.null(minmax)){
    minmax <- range(intmat_trans)
  }
  if (is.null(minmax_for_legend)){
    minmax_for_legend <- range(intmat)
  }
  nr <- build_nr(intmat_trans, minmax = c(minmax[1], minmax[2]), colormap)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      xmin = dt_min, xmax = dt_min,
      ymin = rt_min, ymax = rt_min,
      ggplot2::aes(fill = .data$x),
      data = data.frame(
        x = NA_real_,
        dt_ms_min = dt_min, dt_ms_max = dt_max,
        rt_s_min = rt_min, rt_s_max = rt_max
      )
    ) +
    ggplot2::annotation_raster(
      nr,
      xmin = dt_min, xmax = dt_max,
      ymin = rt_min, ymax = rt_max
    ) +
    scale_fill_gradientn(colors = colormap,   # Custom Spectral colormap
                         limits = minmax_for_legend,              # Use your predefined limits
                         na.value = "#00000000",       # Transparent for NA values
                         trans = trans) +
    ggplot2::lims(
      x = c(dt_min, dt_max),
      y = c(rt_min, rt_max)
    ) +
    ggplot2::labs(
      x = "Drift time (ms)",
      y = "Retention time (s)",
      fill = "Intensity (a.u.)"
    ) +
    ggplot2::theme_minimal()
  print(p)
}

plot_intensity_matrix <- function(intmat, colormap, minmax = NULL, minmax_for_legend = NULL, dt_min = NULL, dt_max = NULL, rt_min = NULL, rt_max = NULL){

  if (is.null(dt_min)) {
    dt_min <- as.numeric(rownames(intmat)[1L])
  }
  if (is.null(dt_max)) {
    dt_max <- as.numeric(rownames(intmat)[nrow(intmat)])
  }
  if (is.null(rt_min)) {
    rt_min <- as.numeric(colnames(intmat)[1L])
  }
  if (is.null(rt_max)) {
    rt_max <- as.numeric(colnames(intmat)[ncol(intmat)])
  }

  trans <- cubic_root_trans()
  intmat_trans <- trans$transform(intmat)
  if (is.null(minmax)){
    minmax <- range(intmat_trans)
  }
  if (is.null(minmax_for_legend)){
    minmax_for_legend <- range(intmat)
  }
  nr <- build_nr(intmat_trans, minmax = c(minmax[1], minmax[2]), colormap)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      xmin = dt_min, xmax = dt_min,
      ymin = rt_min, ymax = rt_min,
      ggplot2::aes(fill = .data$x),
      data = data.frame(
        x = NA_real_,
        dt_ms_min = dt_min, dt_ms_max = dt_max,
        rt_s_min = rt_min, rt_s_max = rt_max
      )
    ) +
    ggplot2::annotation_raster(
      nr,
      xmin = dt_min, xmax = dt_max,
      ymin = rt_min, ymax = rt_max
    ) +
    scale_fill_gradientn(colors = colormap,   # Custom Spectral colormap
                         limits = minmax_for_legend,              # Use your predefined limits
                         na.value = "#00000000",       # Transparent for NA values
                         trans = trans) +
    ggplot2::lims(
      x = c(dt_min, dt_max),
      y = c(rt_min, rt_max)
    ) +
    ggplot2::labs(
      x = "Drift time (ms)",
      y = "Retention time (s)",
      fill = "Intensity (a.u.)"
    ) +
    ggplot2::theme_minimal()
  print(p)
}

# Create colormap for plot_sample()
colormap <- c(
  colorRampPalette(c("black", "#26456E"))(70),          # black → dark blue
  colorRampPalette(c("#26456E", "#1C65A3"))(30),        # dark blue → blue
  colorRampPalette(c("#1C65A3", "#4993C0"))(20),        # blue → light blue
  colorRampPalette(c("#4993C0", "orange"))(20),         # light blue → orange
  colorRampPalette(c("orange", "red"))(60),             # orange → red
  colorRampPalette(c("red", "#CB1618"))(70)             # red → dark red
)

plot_grid <- function(dataset, idx_samples, minmax, num_cols_grid){

  # Crear una lista para almacenar los heatmaps
  heatmaps <- list()
  trans <- cubic_root_trans()

  # Generar los primeros 15 heatmaps y almacenarlos en la lista
  for (i in idx_samples) {
    sample1 <- dataset$getSample(sample = i)
    intmat <- intensity(sample1)
    intmat_trans <- trans$transform(intmat)
    nr <- build_nr(intmat_trans, minmax = c(minmax[1], minmax[2]))
    p <- plot_sample(nr = nr, intmat = intmat, minmax = minmax) + theme(legend.position ="none") + ggtitle(paste0("Sample ", i))
    heatmaps[[i]] <- p
  }

  # Crear un gráfico con una leyenda común usando cowplot
  # Extraer una leyenda de uno de los gráficos
  legend <- get_legend(
    plot_sample(nr = nr, intmat = intmat, minmax = minmax)
  )

  invisible(print(cowplot::plot_grid(
    cowplot::plot_grid(plotlist = heatmaps, ncol = num_cols_grid),
    #legend,
    ncol = 2,
    rel_widths = c(1, 0.1)  # Ajustar el ancho relativo de las columnas
  )))
}

#### Peak list creation ####

create_peak_list <- function(sample_id, object, ...){
  # get sample i form object
  s <- object$getSample(sample_id)

  #### 1) Find peaks in the RIC

  # get the RIC and put it as a GCIMSChromatogram object (so that we can use the findPeaks() later)
  ric <- GCIMSChromatogram(retention_time = rtime(s),
                           intensity = getRIC(s),
                           description = s@description)
  # add peaks information to the GCIMSChromatogram object
  ric <- findPeaks(ric,
                   length_in_xunits = 3,
                   peakwidth_range_xunits = c(10, 50),
                   peakDetectionCWTParams = list(exclude0scaleAmpThresh = TRUE,
                                                 SNR.Th=1, excludeBoundariesSize = 0, ...))
  # save the peaks dataframe in a separate variable
  peaks_rt <- peaks(ric)
  peaks_rt$PeakID <- as.numeric(peaks_rt$PeakID)

  # add 3 empty columns to later save the minimum, the apex and the maximum of the Rt peaks
  peaks_rt <- peaks_rt |>
    dplyr::mutate(dt_peak_min = NA, dt_apex = NA, dt_peak_max = NA)

  #### 2) Find peaks in the TIS

  # get the TIS and put it as a GCIMSSpectrum object (so that we can use the findPeaks() later)
  tis <- getSpectrum(s)
  # intensity 1/3 transformation
  tis@intensity <- sign(tis@intensity)*abs(tis@intensity)^(1/3)
  # add peaks information to the GCIMSSpectrum object
  tis <- findPeaks(tis,
                   length_in_xunits = 0.14,
                   peakwidth_range_xunits = c(0.1, 0.3),
                   peakDetectionCWTParams = list(exclude0scaleAmpThresh = TRUE,
                                                 SNR.Th=1))
  # get the RIP maximum point / get the point where the RIP ends (the RIP is the peak with a higher value of "int_apex_au")
  # "int_apex_au" -> intesnisty apex arbitrary unit
  rip_max <- peaks(tis)$max[which.max(peaks(tis)$int_apex_au)]

  #### 3) Iterate over all the peaks found in Rt

  for (j in 1:nrow(peaks_rt)){
    # for each peak in Rt, get the spectrum in that Rt where the peak is found
    spec_at_specific_rt <- getSpectrum(s, rt_range = peaks_rt$apex[j])
    # intensity 1/3 transformation
    spec_at_specific_rt@intensity <- sign(spec_at_specific_rt@intensity)*abs(spec_at_specific_rt@intensity)^(1/3)
    # add peaks information to the GCIMSSpectrum object
    spec_at_specific_rt <- findPeaks(spec_at_specific_rt,
                                     length_in_xunits = 0.14,
                                     peakwidth_range_xunits = c(0.1, 0.3),
                                     peakDetectionCWTParams = list(exclude0scaleAmpThresh = TRUE, SNR.Th=1))
    # save the peaks dataframe in a separate variable
    peaks_spec_at_specific_rt <- peaks(spec_at_specific_rt)
    # remove the peaks that are before the RIP and the RIP itself
    peaks_spec_at_specific_rt <- peaks_spec_at_specific_rt |>
      dplyr::filter(apex  > rip_max)
    # if there are peaks in Dt at that Rt
    if (nrow(peaks_spec_at_specific_rt) != 0){
      # find the highest peak in the spectrum at the Rt being evaluated
      highest_dt_peak <- which.max(peaks_spec_at_specific_rt$int_apex_au)
      peaks_rt$dt_peak_min[j] <- peaks_spec_at_specific_rt$min[highest_dt_peak]
      peaks_rt$dt_apex[j] <- peaks_spec_at_specific_rt$apex[highest_dt_peak]
      peaks_rt$dt_peak_max[j] <- peaks_spec_at_specific_rt$max[highest_dt_peak]
    }
  }
  # remove the Rt peaks that do not have a peak in Dt
  peaks_rt <- peaks_rt |>
    dplyr::filter(!is.na(dt_apex))
  # save the peaks_rt dataframe in the peak list
  return(peaks_rt)
}

### Find the correspondace between each peak in a sample and the peaks in the reference sample ###

match_peaks <- function(sample_id, peak_list, reference, percentage_movement){
  # get the reference peaks
  reference_peaks <- peak_list[[reference]]
  # get the sample peaks
  sample_peaks <- peak_list[[sample_id]]
  anchors <- data.frame()
  # iterate over all the peaks in each sample
  for (j in 1:nrow(sample_peaks)){
    # find to which peak it corresponds on the reference (it may be more than 1)
    apex_peak_j <- sample_peaks$dt_apex[j]
    peak_match <- which(sample_peaks$dt_apex[j] >= reference_peaks$dt_peak_min &
                          sample_peaks$dt_apex[j] <= reference_peaks$dt_peak_max)
    # if there is a correspondance (either 1 or more peaks of the reference coincide with the peak j)
    if (length(peak_match) > 0){

      # for each of the peaks that were a "match" between peak j and reference
      for (k in 1:length(peak_match)){
        # set the limits in which we want to find the peaks
        min_ref <- reference_peaks[peak_match[k],]$min * (1 - percentage_movement)
        max_ref <- reference_peaks[peak_match[k],]$max * (1 + percentage_movement)

        # if the peak or peaks found are inside the limits
        if (sample_peaks$apex[j] > min_ref & sample_peaks$apex[j] < max_ref){
          # add this peak information
          anchors <- rbind(anchors,
                           data.frame(
                             reference_PeakID = reference_peaks$PeakID[peak_match[k]],
                             reference_rt_apex = reference_peaks$apex[peak_match[k]],
                             reference_dt_apex = reference_peaks$dt_apex[peak_match[k]],
                             sample_PeakID = sample_peaks$PeakID[j],
                             sample_rt_apex = sample_peaks$apex[j],
                             sample_dt_apex = sample_peaks$dt_apex[j]
                           ))

        }
      }
    }
  }
  return(anchors)
}


clean_peak_list <- function(sample_id, peak_list, object, ip_reference){
  s <- object$getSample(sample_id)
  peaks_df <- peak_list[[sample_id]]

  # add a column with the absolute difference between the apex of the reference peak and the apex of the sample peak
  peaks_df <- peaks_df |>
    dplyr::mutate(difference = abs(reference_rt_apex - sample_rt_apex))
  # in case of repeated reference_PeakID, keep just the peak that has a lowest absolute difference between the reference peak and the sample peak
  peaks_df <- peaks_df |>
    dplyr::group_by(reference_PeakID) |>
    dplyr::filter(difference == min(difference)) |>
    dplyr::ungroup()
  # in case of repeated sample_PeakID, keep just the peak that has a lowest absolute difference between the reference peak and the sample peak
  peaks_df <- peaks_df |>
    dplyr::group_by(sample_PeakID) |>
    dplyr::filter(difference == min(difference)) |>
    dplyr::ungroup() |>
    dplyr::select(-difference)

  ip_sample <- rtime(s)[which.min(getRIC(s))]
  ip_row <- c(0, ip_reference, 0, 0, ip_sample, 0)
  peaks_df <- rbind(ip_row, as.matrix(peaks_df))

  return(peaks_df)
}

peaks_alignment <- function(object, reference, percentage_movement, ...){
  sample_ids <- object$sampleNames

  peak_list <- setNames(
    lapply(sample_ids, create_peak_list, object, ...),
    sample_ids)

  peak_list_with_match <- setNames(
    lapply(sample_ids, match_peaks, peak_list, reference = reference, percentage_movement = percentage_movement),
    sample_ids)

  ip_reference <- rtime(object$getSample(reference))[which.min(getRIC(object$getSample(reference)))]

  final_peak_list <- setNames(
    lapply(sample_ids, clean_peak_list, peak_list_with_match, object, ip_reference),
    sample_ids)

  return(final_peak_list)
}

align_pks_sample <- function(object, peak_list){
  rt <- GCIMS::rtime(object)
  int_mat <- GCIMS::intensity(object)
  rt_p <- signal::interp1(x = peak_list[,"sample_rt_apex"], y = peak_list[,"reference_rt_apex"], rt, method = "spline", extrap=TRUE, "monoH.FC")
  object@data <- t(apply(int_mat, 1, signal::interp1, x = rt_p, xi = rt, extrap = TRUE))
  return(object)
}

align_pks_ds <- function(object, reference, percentage_movement, ...){
  peak_list <- peaks_alignment(object, reference, percentage_movement)

  delayed_op <- DelayedOperation(
    name = "align_pks",
    fun = align_pks_sample,
    params_iter = list(peak_list)
  )
  object$appendDelayedOp(delayed_op)
  object$extract_RIC_and_TIS()
  object$extract_dtime_rtime()
  invisible(object)
}


transform_int_s<- function(object){
  int_mat <- GCIMS::intensity(object)
  object@data <- sign(int_mat)*abs(int_mat)^(1/3)
  return(object)
}
transform_int <- function(object){
  delayed_op <- DelayedOperation(
    name = "transform_int",
    fun = transform_int_s
  )
  object$appendDelayedOp(delayed_op)
  object$extract_RIC_and_TIS()
  object$extract_dtime_rtime()
  invisible(object)
}

anti_transform_int_s<- function(object){
  int_mat <- GCIMS::intensity(object)
  object@data <- sign(int_mat)*abs(int_mat)^(3)
  return(object)
}
anti_transform_int <- function(object){
  delayed_op <- DelayedOperation(
    name = "anti_transform_int",
    fun = anti_transform_int_s
  )
  object$appendDelayedOp(delayed_op)
  object$extract_RIC_and_TIS()
  object$extract_dtime_rtime()
  invisible(object)
}


grid_spectrums <- function(dataset, dt_range, rt_range, png_file_name){
  # Crear una lista para almacenar los gráficos
  plots <- list()

  # Generar los gráficos y almacenarlos en la lista
  for (i in 1:55) {
    sample1 <- dataset$getSample(i)
    p <- plot(getSpectrum(sample1, dt_range = dt_range, rt_range = rt_range)) +
      theme(axis.ticks.x=element_blank(),
            axis.title = element_blank(),
            title = element_blank()
      )
    plots[[i]] <- p
  }

  # Calcular el número de filas y columnas para la cuadrícula
  n_col <- 5  # Número de columnas (ajusta según necesites)
  n_row <- ceiling(length(plots) / n_col)

  # Crear la cuadrícula de gráficos
  grid_plot <- marrangeGrob(grobs = plots, nrow = n_row, ncol = n_col)

  png(png_file_name, width = 4000, height = 4000, res = 300)
  grid.arrange(grobs = plots, nrow = n_row, ncol = n_col)
  dev.off()
}

grid_chromatograms <- function(dataset, dt_range, rt_range, png_file_name){
  # Crear una lista para almacenar los gráficos
  plots <- list()

  # Generar los gráficos y almacenarlos en la lista
  for (i in 1:55) {
    sample1 <- dataset$getSample(i)
    p <- plot(getChromatogram(sample1, dt_range = dt_range, rt_range = rt_range)) +
      theme(axis.ticks.x=element_blank(),
            axis.title = element_blank(),
            title = element_blank()
      )
    plots[[i]] <- p
  }

  # Calcular el número de filas y columnas para la cuadrícula
  n_col <- 5  # Número de columnas (ajusta según necesites)
  n_row <- ceiling(length(plots) / n_col)

  # Crear la cuadrícula de gráficos
  grid_plot <- marrangeGrob(grobs = plots, nrow = n_row, ncol = n_col)

  png(png_file_name, width = 4000, height = 4000, res = 300)
  grid.arrange(grobs = plots, nrow = n_row, ncol = n_col)
  dev.off()
}

################################# PLS-DA

folds_plsda_workflow <- function(fold_number, external_folds, X, Y){

  # TRAIN and TEST datasets
  idx_cal <- external_folds[[fold_number]]
  idx_test <- setdiff(1:length(Y), idx_cal)

  Xc <- X[idx_cal,]
  Yc <- Y[idx_cal]

  Xt <- X[idx_test,]
  Yt <- Y[idx_test]

  #### INTERNAL VALIDATION

  # Train model with selected features
  train_model <- mixOmics::splsda(Xc, Yc, ncomp = min(ncol(X), 10))

  perf.splsda.srbct <- mixOmics::perf(train_model, validation = "Mfold",
                                      folds = 5, nrepeat = 20, # use repeated cross-validation
                                      progressBar = FALSE, auc = TRUE) # include AUC values

  # build final train model
  LV <- perf.splsda.srbct$choice.ncomp[5]
  AUC <- perf.splsda.srbct$auc[[LV]][1:3]

  train_final_model <- mixOmics::splsda(Xc, Yc, ncomp = LV)
  # get the AUC of the cross validated train model


  # add the values to the dataframe
  results_hp <- data.frame(ncomp=LV, AUC=I(list(AUC)),
                           ids_train=I(list(idx_cal)), ids_test=I(list(idx_test)))

  ### EXTERNAL VALIDATION

  # Predict
  preds_test <- predict(train_final_model, as.matrix(Xt), dist = "mahalanobis.dist")
  accuracy <- sum(preds_test$class$mahalanobis.dist[,LV] == Yt) / length(Yt)


  # add values to df
  results_hp <- cbind(results_hp, accuracy_external=accuracy, fold=paste0("Fold", fold_number))

  return(results_hp)
}


get_best_num_variables <- function(results_rfe, analysis_name, variable){
  analysis_names <- unique(results_rfe[[analysis_name]])

  name_best_model <- results_rfe %>%
    group_by(!!sym(analysis_name)) %>%
    summarize(median_var = median(!!sym(variable))) %>%
    filter(median_var == min(median_var))
  colnames(name_best_model) <- c("subset", "median")
  name_best_model <- name_best_model %>% pull(subset)
  return(name_best_model)
}


select_variables <- function(results_rfe,analysis_name, name_best_model){
  models_list <- results_rfe %>% filter(!!sym(analysis_name)==name_best_model)  %>% pull(variables)
  sorted_table <- sort(table(unlist(models_list)), decreasing=TRUE)
  top_variables <- names(sorted_table)[1:as.integer(name_best_model)]
  return(top_variables)
}

recursive_feature_elimination_PLSDA <- function(rep, Xc, Yc){

  # build the dataframe to store the values
  results_hp <- data.frame()

  # initialization
  num_vars_to_keep <- ncol(Xc)
  vars_to_keep <- colnames(Xc)

  while (num_vars_to_keep >= 1){
    print(num_vars_to_keep)

    # select how many components we want to test
    ncomp <- min(num_vars_to_keep, 5)
    # build the PLSDA object
    plsda_obj <- mixOmics::splsda(Xc[,vars_to_keep], Yc, ncomp = ncomp)

    # undergo performance evaluation in order to tune the number of components to use
    perf.splsda.srbct <- mixOmics::perf(plsda_obj, validation = "Mfold",
                                        folds = 5, nrepeat = 5, # use repeated cross-validation
                                        progressBar = FALSE, auc = TRUE) # include AUC values

    # # select the optimal number of components
    # optimal_ncomp <- perf.splsda.srbct$choice.ncomp[1]
    # # build final model
    # plsda_obj_final_model <- mixOmics::splsda(X[,vars_to_keep], Y, ncomp = optimal_ncomp)
    #
    # # get the AUC of the cross validated final model
    # auc <- perf.splsda.srbct$auc[optimal_ncomp][[1]][["AUC.mean"]]
    LV <- perf.splsda.srbct$choice.ncomp[5]
    AUC <- perf.splsda.srbct$auc[[LV]][1:3]
    error_rate <- perf.splsda.srbct[["error.rate"]]$overall[LV,3]

    plsda_obj_final_model <- mixOmics::splsda(Xc[,vars_to_keep], Yc, ncomp = LV)

    # add the values to the dataframe
    results_hp <- rbind(results_hp, data.frame(ncomp=LV, error_rate = error_rate, AUC=I(list(AUC)),
                                               num_vars=length(vars_to_keep), variables=I(list(vars_to_keep)), rep=rep))

    # select the variables to keep for the next step of the RFE
    num_vars_to_keep <- round(num_vars_to_keep*80/100) # number of variables to keep
    vips_final_model <- names(sort(mixOmics::vip(plsda_obj_final_model)[,LV], decreasing = TRUE))
    vars_to_keep <- vips_final_model[1:num_vars_to_keep]

    if (num_vars_to_keep <= 2) {
      num_vars_to_keep <- num_vars_to_keep - 1

    }
  }
  return(results_hp)
}


folds_plsda_workflow_RFE <- function(fold_number, external_folds, X, Y){

  # TRAIN and TEST datasets
  idx_cal <- external_folds[[fold_number]]
  idx_test <- setdiff(1:nrow(X), idx_cal)

  Xc <- X[idx_cal,]
  Yc <- Y[idx_cal]

  Xt <- X[idx_test,]
  Yt <- Y[idx_test]

  ### RFE

  print(paste0("Starting RFE from Fold = ", fold_number))
  rfe_res <- recursive_feature_elimination_PLSDA(1, Xc, Yc)

  # Feature selection
  best_num_vars <- get_best_num_variables(rfe_res,"num_vars", "error_rate")[1]
  vars_to_keep <- select_variables(rfe_res,"num_vars", best_num_vars)
  LV <- rfe_res %>% filter(num_vars==best_num_vars) %>% pull(ncomp)
  error_rate_train <- rfe_res %>% filter(num_vars==best_num_vars) %>% pull(error_rate)
  train_final_model <- mixOmics::splsda(Xc[,vars_to_keep], Yc, ncomp = LV)

  ### EXTERNAL VALIDATION

  # Predict
  preds_test <- predict(train_final_model, as.matrix(Xt[,vars_to_keep]), dist = "mahalanobis.dist")
  accuracy <- sum(preds_test$class$mahalanobis.dist[,LV] == Yt) / length(Yt)


  # add values to df
  results_hp <- data.frame(fold=paste0("Fold", fold_number), error_rate_train = error_rate_train, accuracy_external=accuracy,
                           num_vars = best_num_vars, LV = LV, variables = I(list(vars_to_keep)), ids_train=I(list(idx_cal)), ids_test=I(list(idx_test)))

  return(results_hp)
}

folds_plsda_workflow_vip1 <- function(fold_number, external_folds, X, Y){

  # TRAIN and TEST datasets
  idx_cal <- external_folds[[fold_number]]
  idx_test <- setdiff(1:length(Y), idx_cal)

  Xc <- X[idx_cal,]
  Yc <- Y[idx_cal]

  Xt <- X[idx_test,]
  Yt <- Y[idx_test]

  #### INTERNAL VALIDATION

  # Train model with selected features
  train_model <- mixOmics::splsda(Xc, Yc, ncomp = min(ncol(X), 10))

  perf.splsda.srbct <- mixOmics::perf(train_model, validation = "Mfold",
                                      folds = 5, nrepeat = 20, # use repeated cross-validation
                                      progressBar = FALSE, auc = TRUE) # include AUC values

  # build final train model
  LV <- perf.splsda.srbct$choice.ncomp[5]
  AUC <- perf.splsda.srbct$auc[[LV]][1:3]

  train_final_model <- mixOmics::splsda(Xc, Yc, ncomp = LV)

  vips_final_model <- mixOmics::vip(train_final_model)[,LV]
  vars_vips_greater_1 <- names(subset(vips_final_model, vips_final_model > 1))

  train_final_model <- mixOmics::splsda(Xc[, vars_vips_greater_1], Yc, ncomp = LV)

  # add the values to the dataframe
  results_hp <- data.frame(ncomp=LV, AUC=I(list(AUC)),
                           ids_train=I(list(idx_cal)), ids_test=I(list(idx_test)), variables=I(list(vars_vips_greater_1)))

  ### EXTERNAL VALIDATION
  # Predict
  preds_test <- predict(train_final_model, as.matrix(Xt[, vars_vips_greater_1]), dist = "mahalanobis.dist")
  accuracy <- sum(preds_test$class$mahalanobis.dist[,LV] == Yt) / length(Yt)


  # add values to df
  results_hp <- cbind(results_hp, accuracy_external=accuracy, fold=paste0("Fold", fold_number))

  return(results_hp)
}













