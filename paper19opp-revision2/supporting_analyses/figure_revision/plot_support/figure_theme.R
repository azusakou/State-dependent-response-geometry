palette_control <- c(
  neutral_dark = "#2B2B2B",
  neutral_mid = "#7A7A7A",
  neutral_light = "#D9D9D9",
  genetic = "#4C78A8",
  combination = "#7B6BA8",
  multiome = "#2A9D8F",
  dynamics = "#E3A857",
  pass = "#2A9D8F",
  pending = "#A0A0A0",
  fail = "#C95A4B"
)

theme_nature_control <- function(base_size = 7, base_family = "Helvetica") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.4),
      legend.title = ggplot2::element_text(size = base_size - 0.2),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      strip.text = ggplot2::element_text(size = base_size, face = "plain"),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(7, 8, 7, 7),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 95,
                       dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  # Enforce the no-panel-title contract even when a plot overrides the theme.
  if (inherits(plot, "patchwork")) {
    plot <- plot & ggplot2::theme(plot.title = ggplot2::element_blank(),
                                  plot.subtitle = ggplot2::element_blank())
  } else {
    plot <- plot + ggplot2::theme(plot.title = ggplot2::element_blank(),
                                  plot.subtitle = ggplot2::element_blank())
  }

  svglite::svglite(paste0(filename, ".svg"), width = width_in,
                   height = height_in)
  print(plot)
  grDevices::dev.off()

  # Quartz embeds TrueType fonts on macOS without an XQuartz dependency.
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    grDevices::quartz(type = "pdf", file = paste0(filename, ".pdf"),
                     width = width_in, height = height_in, family = "Helvetica")
  } else {
    grDevices::cairo_pdf(paste0(filename, ".pdf"), width = width_in,
                         height = height_in, family = "Helvetica")
  }
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(paste0(filename, ".tiff"), width = width_in,
                 height = height_in, units = "in", res = dpi,
                 compression = "lzw")
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(paste0(filename, "_preview.png"), width = width_in,
                height = height_in, units = "in", res = 300)
  print(plot)
  grDevices::dev.off()
}
