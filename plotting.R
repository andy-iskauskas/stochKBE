####################
# Plotting Helpers #
####################

library(ggplot2)
library(viridis)

## Labellers for prettifying plots
plot_labels <- as_labeller(
  c(
    "Vbulk" = "Bulk Data Only",
    "Vbound" = "Boundary Data Only",
    "Vboth" = "Bulk and Boundary Data",
    "Vno" = "Prior (no data)",
    "Ebulk" = "Bulk Data Only",
    "Ebound" = "Boundary Data Only",
    "Eboth" = "Bulk and Boundary Data",
    "Eno" = "Prior (no data)",
    "Ibulk" = "Bulk Data Only",
    "Ibound" = "Boundary Data Only",
    "Iboth" = "Bulk and Boundary Data",
    "Ino" = "Prior (no data)"
  )
)

## Implausibility contour definition
imp_breaks <- c(0, 0.3, 0.7, 1, 1.3, 1.7, 2, 2.3, 2.7,
                3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 15, Inf)
imp_names <- c(0, '', '', 1, '', '', 2, '', '', 3, '',
               '', '', 5, '', '', '', 10, 15, '')
redgreen <- c('#00FF00', '#18FF00', '#31FF00', '#49FF00', '#62FF00',
              '#7AFF00', '#93FF00', '#ABFF00', '#C4FF00', '#DDFF00',
              '#E0E200', '#E4C600', '#E8AA00', '#EC8D00', '#EF7100',
              '#F35500', '#F73800', '#FB1C00', '#FF0000', '#FF0000')

## Plot function for shorthand plotting
# Takes data to plot, data prefix to identify columns, the omega value to slice at,
# a title prefix, training points (if wanted to overlay those), and a distance from
# the slice beyond which points are plotted in grey
grid_plot <- function(data, prefix, omega_pt, title_add = "", training_pts = NULL, pt_dist = Inf) {
  if (prefix == "E") p_title <- paste(title_add, "Emulator Expectation")
  if (prefix == "V") p_title <- paste(title_add, "Emulator Variance")
  if (prefix == "I") p_title <- paste(title_add, "Emulator Implausibility")
  g <- ggplot(data = subset(data[grep(prefix, data$name),], omega == omega_pt),
              aes(x = beta, y = gamma))
  if (prefix != "I") {
    g <- g + geom_contour_filled(aes(z = value)) +
      scale_fill_viridis(discrete = TRUE) + guides(fill = guide_legend(ncol = 1))
  }
  else {
    g <- g + geom_contour_filled(aes(z = value), colour = 'black', linewidth = 0.1, breaks = imp_breaks) +
      geom_contour(aes(z = value), breaks = c(0, 3, Inf), colour = 'black') + 
      scale_fill_manual(values = redgreen, name = "I", labels = imp_names,
                               guide = guide_legend(ncol = 1, reverse = TRUE))
  }
  if (!is.null(training_pts)) {
    g <- g + geom_point(data = subset(training_pts, abs(omega - omega_pt) < pt_dist)) +
      geom_point(data = subset(training_pts, abs(omega - omega_pt) > pt_dist), col = "grey40",
                 pch = 4, cex = 0.5)
  }
  g <- g + facet_wrap(vars(name), nrow = 2, ncol = 2, labeller = plot_labels) +
    labs(title = p_title)
  return(g)
}