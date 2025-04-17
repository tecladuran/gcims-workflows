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
