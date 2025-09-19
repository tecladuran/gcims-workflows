library(patchwork)
library(cowplot)

plot_ric_2 <- function(urine) {
  ric <- getRIC(urine)
  melt <- reshape2::melt(ric)
  custom_colors <- paletteer::paletteer_c("ggthemes::Blue-Teal", nrow(ric))  # nombre de mostres

  ggplot() +
    geom_line(
      data = melt,
      aes(x = retention_time_s, y = value, group = SampleID, color = SampleID),
      linewidth = 0.27
    ) +
    scale_color_manual(values = custom_colors) +
    labs(
      x = "Retention time (s)",
      y = "RIC Intensity (a.u.)",
      color = "SampleID"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}


plot_pca <- function(data, color_var, title, color_scale) {
  ggplot(data, aes_string(x = "PC1", y = "PC2", color = color_var)) +
    geom_point(size = 2) +
    color_scale +
    theme_minimal(base_size = 16)+
    labs(
      title = title,
      x = axis_labels[1],
      y = axis_labels[2],
      color = gsub("_", " ", color_var)   # Elimina els underscores
    ) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5)
    )
}

plot_pca_2 <- function(data1, data2, color_var, title1, title2, color_scale) {

  xlim <- range(data1$PC1, data2$PC1)
  ylim <- range(data1$PC2, data2$PC2)

  p1 <- ggplot(data1, aes_string(x = "PC1", y = "PC2", color = color_var)) +
    geom_point(size = 2) +
    color_scale +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_minimal(base_size = 13) +
    labs(title = title1, x = axis_labels[1], y = axis_labels[2]) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "right"
    )

  p2 <- ggplot(data2, aes_string(x = "PC1", y = "PC2", color = color_var)) +
    geom_point(size = 2) +
    color_scale +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_minimal(base_size = 13) +
    labs(title = title2, x = axis_labels[1], y = axis_labels[2]) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      legend.position = "right"
    )

  return((p1 + p2) + plot_layout(guides = "collect"))
}
