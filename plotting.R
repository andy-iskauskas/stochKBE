####################
# Plotting Helpers #
####################

library(ggplot2)
library(viridis)
library(latex2exp)

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
grid_plot <- function(data, prefix, plot_names = c("beta", "gamma"), title_add = "", training_pts = NULL,
                      breaks = NULL, labels = NULL, viridoption = "A") {
  data$name <- factor(data$name, levels = paste0(prefix, c("no", "bound", "bulk", "both")))
  if (prefix == "E") p_title <- paste(title_add, "Emulator Expectation")
  if (prefix == "V") p_title <- paste(title_add, "Emulator Variance")
  if (prefix == "I") p_title <- paste(title_add, "Emulator Implausibility")
  dat_subs <- data[grep(prefix, data$name),]
  g <- ggplot(data = dat_subs, aes(x = .data[[plot_names[1]]], y = .data[[plot_names[2]]]))
  if (prefix != "I") {
    if (is.null(breaks)) {
      g <- g + geom_contour_filled(aes(z = value)) +
        scale_fill_viridis(discrete = TRUE, option = viridoption) + guides(fill = guide_legend(ncol = 1))
    }
    else {
      if (is.null(labels))
        g <- g + geom_contour_filled(aes(z = value), breaks = breaks) +
          scale_fill_viridis(discrete = TRUE, option = viridoption) + guides(fill = guide_legend(ncol = 1))
      else
        g <- g + geom_contour_filled(aes(z = value), breaks = breaks) +
          scale_fill_viridis(discrete = TRUE, option = viridoption, labels = labels) + guides(fill = guide_legend(ncol = 1))
    }
  }
  else {
    g <- g + geom_contour_filled(aes(z = value), colour = 'black', linewidth = 0.1, breaks = imp_breaks) +
      geom_contour(aes(z = value), breaks = c(0, 3, Inf), colour = 'black') + 
      scale_fill_manual(values = redgreen, name = "I", labels = imp_names,
                               guide = guide_legend(ncol = 1, reverse = TRUE))
  }
  if (!is.null(training_pts)) {
    g <- g + geom_point(data = training_pts)
  }
  g <- g + facet_wrap(vars(name), nrow = 2, ncol = 2, labeller = plot_labels) +
    labs(title = p_title, x = TeX(paste0("$\\", plot_names[[1]], "$")),
         y = TeX(paste0("$\\", plot_names[[2]], "$")))
  return(g)
}

comp_labeller <- as_labeller(
  c(
    "Old" = "20-point LHD",
    "Naive" = "40-point augmented LHD",
    "Uniform" = "Uniform Repetition IMSPE",
    "New" = "Full IMSPE"
  ))
comparison_plot <- function(data, facet_names, plot_names = c("beta", "gamma"), 
                            levels_name = "level", breaks = NULL, labels = NULL,
                            viridoption = "D") {
  data_reshape <- tidyr::pivot_longer(data, cols = all_of(facet_names))
  data_reshape$name <- factor(data_reshape$name, levels = facet_names)
  g <- ggplot(data = data_reshape, aes(x = .data[[plot_names[1]]], y = .data[[plot_names[2]]], z = value))
  if (!is.null(breaks)) {
    g <- g + geom_contour_filled(breaks = breaks)
    if (!is.null(labels))
      g <- g + scale_fill_viridis(discrete = TRUE, labels = labels, name = levels_name,
                                  option = viridoption)
    else
      g <- g + scale_fill_viridis(discrete = TRUE, name = levels_name,
                                  option = viridoption)
  }
  else
    g <- g + geom_contour_filled() + scale_fill_viridis(discrete = TRUE, name = levels_name,
                                                        option = viridoption)
  return(g + labs(x = TeX(paste0("$\\", plot_names[[1]], "$")), y = TeX(paste0("$\\", plot_names[[2]], "$"))) +
           guides(fill = guide_legend(ncol = 1)) + facet_wrap(vars(name), nrow = 2, labeller = comp_labeller))
}
