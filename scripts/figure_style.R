# Common style for the manuscript figures, R side (same specification as scripts/figure_style.py; user, 2026-11-24).
# width 7.1 in, height <= 9 in, 600 dpi; Arial 7 pt (ticks 6.5, small 6, panel tags 9 bold); axis lines 0.6 pt;
# no box around panels, no minor grid; top / right axis lines off.
# ggplot2 linewidth is in mm: 1 pt = 0.3528 mm.
PT <- 1 / 2.845276
theme_paper <- function(grid_y = FALSE) {
  ggplot2::theme_classic(base_size = 7, base_family = "Arial") +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.6 * PT, colour = "#111111"),
                   axis.ticks = ggplot2::element_line(linewidth = 0.6 * PT, colour = "#111111"),
                   axis.ticks.length = grid::unit(2.5, "pt"),
                   axis.text = ggplot2::element_text(size = 6.5, colour = "#111111"),
                   legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 6.5),
                   plot.tag = ggplot2::element_text(face = "bold", size = 9),
                   panel.grid.major.y = if (grid_y) ggplot2::element_line(linewidth = 0.4 * PT, colour = "#e6e5e1")
                                        else ggplot2::element_blank())
}
