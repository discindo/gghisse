#' @name gghisse
#' @docType package
#' @keywords package
#'
#' @import hisse dplyr ggplot2
#'
#' @importFrom stats median quantile sd
#' @importFrom stringr str_replace str_extract str_detect str_remove regex
#' @importFrom ggtree ggtree geom_tiplab geom_tiplab2
#' @importFrom treeio as.treedata
#' @importFrom ggridges geom_density_ridges
#' @importFrom colorplaner scale_color_colorplane interpolate_projection guide_colorplane
#' @importFrom viridis viridis
#' @importFrom ape branching.times
#' @importFrom purrr map
#' @importFrom rlang .data
#' @importFrom shiny tags 
#'
NULL


##### --- TO DO ----------------------------------- #####

# 5. Support region functions. Making tables, plots, and contour plots?
# 6. g_ functions for HiGeoSSE
# 7. diagram functions

##### --- HiSSE functions ------------------------- #####

#' Prepare HiSSE marginal reconstruction for plotting
#'
#'@description A function to process the output from \code{hisse::MarginRecon}. First we get model averaged rates for tips and nodes of the phylogeny. Then convert the tree from \code{phylo} format into a \code{treedata} class and assign the data.frame of mode-averaged rates to the \code{@data} slot. We output the node and tip rates tables and the reformated tree object for downstream plotting with \code{ggtree} and \code{ggplot}.
#'
#'@param hisse_recon An object produced by \code{hisse::MarginRecon}, or a list of such objects over which we can model-average using \code{hisse::GetModelAveRates}
#'
#'@return A list with three components: \itemize{
#'  \item \code{tree_data} the phylogeny with associated data object with reconstructed ancestral, model-averaged node and tip states and rates
#'  \item \code{tip_rates/node_rates} a tibble of reconstructed ancestral, model-averaged tip/node states and rates.
#'  }
#'
#'
#'@examples
#'
#'data("diatoms")
#'processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)
#'processed_hisse$tip_rates
#'processed_hisse$node_rates
#'processed_hisse$tree_data
#'
#'@export

h_process_recon <- function(hisse_recon) {
  tip.rates <-
    GetModelAveRates(x = hisse_recon, type = "tips") %>%
    as_tibble()
  colnames(tip.rates)[1] <- "id"
  tip.rates$id <- as.character(tip.rates$id)

  nod.rates <-
    GetModelAveRates(x = hisse_recon, type = "nodes") %>%
    as_tibble()

  nod.rates$id <- as.character(nod.rates$id)
  tip.rates$id <- as.character(tip.rates$id)
  both.rates <- bind_rows(tip.rates, nod.rates)

  if (class(hisse_recon) == "list") {
    tree <- as.treedata(hisse_recon[[1]]$phy)
  }

  if (class(hisse_recon) == "hisse.states") {
    tree <- as.treedata(hisse_recon$phy)
  }

  tree@data <- both.rates

  return(list(
    "tree_data" = tree,
    "tip_rates" = tip.rates,
    "node_rates" = nod.rates
  ))
}

#' Plot diversification rates estimated by a HiSSE model (HiSSE, BiSSE, CID2, CID4)
#'
#' @description A function to plot a scatterplot of (model-averaged) diversification rates in the alternative states. We can change the rate plotted on the y axis and modify the label for the x-axis (your binary 0/1 trait).
#'
#' @param processed_recon An object produced by \code{h_process_recon}
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param states_names A character vector of length two giving the translation for states 0 and 1
#' @param colors Colors for the points in the two alternate states
#' @param plot_as_waiting_time Logical. Whether to plot the rates or their inverse (waiting times)
#'
#' @return A jittered scatterplot of (model averaged) tip-associated rates.
#'
#' @examples
#'
#'data("diatoms")
#'processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)
#'hisse_rates_plot <- h_scatterplot(
#'  processed_recon=processed_hisse,
#'  parameter="turnover")
#'
#'# modifications are easy with ggplot
#'# change x axis tick labels
#'hisse_rates_plot +
#'  scale_x_discrete(breaks=c(0,1), labels=c("plankton", "benthos"))
#'
#'# change the position of the legend
#'hisse_rates_plot +
#'  scale_x_discrete(breaks=c(0,1), labels=c("plankton", "benthos")) +
#'  theme(legend.position="top")
#'
#'# use expressions (for greek letters) in the axis labels
#'hisse_rates_plot +
#'  scale_x_discrete(breaks=c(0,1), labels=c("plankton", "benthos")) +
#'  theme(legend.position="top") +
#'  labs(y=expression(paste(tau, "=", lambda, "+", mu))) +
#'  theme(axis.text.y=element_text(size=15))
#'
#'@export

h_scatterplot <-
  function(processed_recon,
           parameter = "turnover",
           states_names = c("0", "1"),
           colors = c("yellow", "violet"),
           plot_as_waiting_time = FALSE) {
    tip.rates <- processed_recon$tip_rates
    tip.rates <- tip.rates %>%
      mutate("f_state" = ifelse(.data$state == 0, states_names[1], states_names[2])) %>%
      mutate("f_state" = factor(.data$f_state))

    if (plot_as_waiting_time) {
      tip.rates <- mutate(tip.rates, wanted = 1 / !!as.name(parameter))
    } else {
      tip.rates <- mutate(tip.rates, wanted = !!as.name(parameter))
    }

    tip.rates.sum <- tip.rates %>%
      group_by(.data$f_state) %>%
      select(.data$f_state, .data$wanted) %>%
      summarise_if(.predicate = is.numeric,
                   .funs = list(
                     "Mean" = mean,
                     "SD" = sd,
                     "Max" = max
                   )) %>%
      mutate("wanted" = .data$Mean) #to bypass warning in geom_errorbar

    result <-
      ggplot(data = tip.rates,
             aes(
               x = .data$f_state,
               y = .data$wanted,
               # color = .data$f_state,
               fill = .data$f_state
             )) +
      geom_point(
        alpha = .7,
        size = 2,
        position = position_jitter(width = .15, height = 0),
        pch=21
      ) +
      geom_errorbar(
        data = tip.rates.sum,
        width = .05,
        aes(
          x = .data$f_state,
          ymax = .data$Mean + .data$SD,
          ymin = .data$Mean - .data$SD
        ),
        position = position_nudge(x = 0.3, y = 0)
      ) +
      geom_point(
        data = tip.rates.sum,
        size = 2.5,
        position = position_nudge(x = 0.3, y = 0),
        aes(x = .data$f_state, y = .data$Mean),
        pch = 21
      ) +
      scale_color_manual(values = colors, name = "") +
      scale_fill_manual(values = colors, name = "") +
      scale_y_continuous(breaks = pretty(c(0, max(
        tip.rates.sum$Max
      )), n = 8)) +
      theme_classic() +
      theme(legend.position = "right",
            legend.key.size = unit(x = .6, units = "cm")) +
      labs(x = "", y = parameter) +
      theme_xy
    return(result)
  }

#' Plot diversification rates estimated by a HiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a dotplot of (model-averaged) diversification rates in the alternative states.
#'
#' @param processed_recon An object produced by \code{h_process_recon}
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param states_names A character vector of length two giving the translation for states 0 and 1.
#' @param bin_width The width of bins for the dotplot. Treat this as any histogram. Testing several different bin width values is recommended.
#' @param colors Colors for the points in the two alternate states
#' @param plot_as_waiting_time Whether to plot the rates or their inverse (waiting times)
#'
#' @return A dotplot of tip-associated rates (possibly model averaged).
#'
#' @examples
#'
#'data("diatoms")
#'processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)
#'
#' h_dotplot(
#'   processed_recon = processed_hisse,
#'   parameter = "turnover",
#'   states_names = c("Plankton", "Benthos"),
#'   bin_width = 0.2,
#'   colors=c("orange", "violet"),
#'   plot_as_waiting_time = TRUE
#' ) + labs(x = "", y = "Waiting time (My)", title = "Turnover")
#'
#'# see ?h_scatterplot for examples for modifying the graph using ggplot2
#'
#' @export

h_dotplot <-
  function(processed_recon,
           parameter = "turnover",
           states_names = c("Marine", "Freshwater"),
           bin_width = 0.1,
           colors,
           plot_as_waiting_time = TRUE) {
    tip.rates <- processed_recon$tip_rates
    tip.rates <- tip.rates %>%
      mutate("f_state" = ifelse(.data$state == 0, states_names[1], states_names[2])) %>%
      mutate("f_state" = factor(.data$f_state))

    if (plot_as_waiting_time) {
      tip.rates <- mutate(tip.rates, "wanted" = 1 / !!as.name(parameter))
    } else {
      tip.rates <- mutate(tip.rates, "wanted" = !!as.name(parameter))
    }

    tip.rates.sum <- tip.rates %>%
      group_by(.data$f_state) %>%
      select(.data$f_state, .data$wanted) %>%
      summarise_if(.predicate = is.numeric,
                   .funs = list(
                     "Mean" = mean,
                     "SD" = sd,
                     "Max" = max
                   )) %>%
      mutate("wanted" = .data$Mean) #to bypass warning in geom_errorbar

    sss <-
      ggplot(data = tip.rates,
             aes(
               x = .data$f_state,
               y = .data$wanted,
               fill = .data$f_state
             )) +
      geom_dotplot(
        alpha = .75,
        colour = "white",
        binwidth = bin_width,
        dotsize = 1.2,
        stackratio = 0.5,
        binpositions = "all",
        binaxis = "y",
        stackdir = "up"
      ) +
      geom_errorbar(
        data = tip.rates.sum,
        width = .05,
        aes(
          x = .data$f_state,
          ymax = .data$Mean + .data$SD,
          ymin = .data$Mean - .data$SD
        ),
        position = position_nudge(x = -0.1, y = 0)
      ) +
      geom_point(
        data = tip.rates.sum,
        size = 2.5,
        position = position_nudge(x = -0.1, y = 0),
        aes(x = .data$f_state, y = .data$Mean),
        pch = 21
      ) +
      scale_color_manual(values = colors, name = "") +
      scale_fill_manual(values = colors, name = "") +
      scale_y_continuous(breaks = pretty(c(0, max(
        tip.rates.sum$Max
      )), n = 8)) +
      theme_classic() +
      theme(
        legend.position = "none",
        legend.background = element_blank(),
        legend.key.size = unit(x = .6, units = "cm")
      ) +
      labs(x = "", y = parameter)+
      theme_xy
    return(sss)
  }

#' Plot diversification rates estimated by a HiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a ridgeline of (model-averaged) diversification rates in the alternative states. We can change the rate plotted on the y axis and modify the label for the x-axis (your binary 0/1 trait). Further modifications are straightforward with ggplot (see examples).
#'
#' @param processed_recon An object produced by \code{h_process_recon}
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param states_names The names for character states
#' @param fill_colors Colors for the density polygons
#' @param line_colors Colors for the lines and points
#' @param plot_as_waiting_time Whether to plot the rates (FALSE, default) or their inverse (waiting times)
#'
#' @return A ridgeline plot of tip-associated rates (possibly model averaged).
#'
#' @examples
#'
#'data("diatoms")
#'processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)
#'
#'h_ridgelines(
#'  processed_recon = processed_hisse,
#'  states_names = c("Plankton", "Benthos"),
#'  parameter = "extinction",
#'  line_colors = c("black", "black"),
#'  fill_colors = c("orange", "violet"))
#'
#'@export

h_ridgelines <- function(processed_recon,
                         parameter = "turnover",
                         states_names = c("Marine", "Freshwater"),
                         line_colors = c("black", "black"),
                         fill_colors = c("yellow", "red"),
                         plot_as_waiting_time = FALSE) {
  tip.rates <- processed_recon$tip_rates
  tip.rates$f_state <-
    factor(ifelse(tip.rates$state == 0,
                  states_names[1],
                  states_names[2]),
           levels = states_names)

  if (plot_as_waiting_time) {
    tip.rates <- mutate(tip.rates, "wanted" = 1 / !!as.name(parameter))
  } else {
    tip.rates <- mutate(tip.rates, "wanted" = !!as.name(parameter))
  }

  tip.rates.sum <- tip.rates %>%
    group_by(.data$f_state) %>%
    select(.data$f_state, .data$wanted) %>%
    summarise_at(
      .vars = vars(.data$wanted),
      .funs = funs(
        "Mean" = mean,
        "SD" = sd,
        "Median" = median,
        # "q95" = quantile(., 0.95),
        # "q05" = quantile(., 0.05),
        "Min" = min,
        "Max" = max
      )
    )
  print(tip.rates.sum)

  ggplot(data = tip.rates,
         aes(
           x = .data$wanted,
           y = .data$f_state,
           fill = .data$f_state,
           color = .data$f_state
         )) +
    geom_density_ridges(
      alpha = 0.75,
      size = 0.75,
      jittered_points = TRUE,
      scale = .95,
      rel_min_height = .01,
      point_shape = "|",
      point_size = 1,
      position = position_nudge(y = rep(-.2, 2))
    ) +
    geom_errorbarh(
      data = tip.rates.sum,
      position = position_nudge(y = rep(-0.3, 2)),
      aes(
        xmin = .data$Mean - .data$SD,
        xmax = .data$Mean + .data$SD,
        y = .data$f_state,
        colour = .data$f_state
      ),
      height = 0.05,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = tip.rates.sum,
      pch = 21,
      position = position_nudge(y = rep(-0.3, 2)),
      aes(
        y = .data$f_state,
        x = .data$Mean,
        colour = .data$f_state
      ),
      size = 3,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(name="", values = fill_colors) +
    scale_color_manual(name="", values = line_colors) +
    scale_x_continuous(breaks = pretty(c(0, max(
      tip.rates.sum$Max
    )), n = 8)) +
    labs(y = "", x = parameter) +
    theme_classic() +
    theme(legend.position = "none")+
    theme_xy
}

#' Plot HiSSE model-averaged marginal ancestral state reconstruction for the trait
#'
#' @description A function to plot a (model-averaged) marginal ancestral reconstruction for the trait data.
#'
#' @param processed_recon An object produced by \code{h_process_recon}
#' @param show_tip_labels Logical, whether to plot tip labels. Default is FALSE because it is difficult to plot legible tip labels for larger trees common in this type of analysis. See \code{?m_trait_recon} for a good manual solution.
#' @param trait_name The name of the trait to be used for guide title
#' @param states_names The names for character states
#' @param colors Colors for the two character states. If \code{discrete=FALSE} these colors will be used as the min/max of a color gradient.
#' @param discrete Logical. Whether to discretize the probabilities of ancestral states into binary (0/1)
#' @param cutoff A decimal to be used as a threshold for discretizing
#' @param tree_layout A layout for the tree. Available options are 'rectangular' (default), 'slanted', 'circular', 'fan' and 'radial'.
#' @param tree_direction 'right' (default), 'left', 'up', or 'down' for rectangular and slanted tree layouts
#' @param time_axis_ticks numeric giving the number of ticks for the time axis (default=10)
#' @param open_angle The degrees of empty space between the first and last tip. Only works for \code{tree_layout = 'fan'} and allows for a little more space around axis tick labels.
#'
#' @return A plot of the phylogeny with branches colored by hisse-inferred marginal ancestral states.
#'
#' @examples
#'
#'data("diatoms")
#'processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)
#'
#'map_continuous <-
#'  h_trait_recon(
#'    processed_recon = processed_hisse,
#'    trait_name = "", discrete=FALSE, cutoff=.5)
#'
#'# change colors, your can pass the trait name to `name=` to title the colorbar
#'map_continuous + scale_color_gradient(name="", low = "#132B43", high = "#56B1F7")
#'
#'map_discrete <-
#'  h_trait_recon(
#'    processed_recon = processed_hisse,
#'    trait_name = "", discrete=TRUE, cutoff=.5)
#'
#'# change colors
#'map_discrete + scale_color_manual(name="", values = c("red", "blue"))
#'
#'@export

h_trait_recon <-
  function(processed_recon,
           show_tip_labels = FALSE,
           trait_name = "trait",
           states_names = c("0", "1"),
           colors= c("orange", "violet"),
           discrete = FALSE,
           cutoff = 0.5,
           tree_layout = "rectangular",
           tree_direction = "right",
           time_axis_ticks = 10,
           open_angle = 10) {
    if (!tree_layout %in% c('rectangular', 'circular', 'slanted', 'fan', 'radial')) {
      stop("The selected tree layout is not supported.")
    }

    tree <- processed_recon$tree_data@phylo
    state <- processed_recon$tree_data@data$state
    agemax <- tree %>% branching.times() %>% max()

    ggg <-
      ggtree(
        tr = tree,
        layout = tree_layout,
        size = 0.45,
        open.angle = open_angle
      ) +
      theme(
        # legend.direction = "horizontal",
        legend.position = "right",
        legend.key.size = unit(x = .5, units = "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(),
        plot.background = element_blank()
      )

    if (discrete) {
      d_state <- ifelse(state > cutoff, 1, 0)
      d_state <- ifelse(d_state == 0, states_names[1], states_names[2]) %>% as.factor()
      ggg <- ggg +
        aes(color = d_state) +
        scale_color_manual(name = trait_name, values = colors) +
        guides(color = guide_legend(override.aes = list(size = 4)))
    } else {
      d_state <- state
      ggg <- ggg +
        aes(color = d_state) +
        scale_color_gradient(name = trait_name, low = colors[1], high = colors[2])
    }

    ggg <-
      tree_flip(
        ggtree_object = ggg,
        show_tip_labels = show_tip_labels,
        tree_layout = tree_layout,
        tree_direction = tree_direction,
        time_axis_ticks = time_axis_ticks,
        agemax = agemax
      )

    return(ggg + theme(plot.margin = unit(rep(.1, 4), "in")))
  }

#' Plot HiSSE model-averaged marginal ancestral state reconstruction for diversification rates
#'
#' @description A function to plot a (model-averaged) marginal ancestral reconstruction for the estimated diversification rates.
#'
#' @param processed_recon An object produced by \code{h_process_recon}
#' @param show_tip_labels Logical, whether to plot tip labels. Default is FALSE because it is difficult to plot legible tip labels for larger trees common in this type of analysis. See \code{?m_trait_recon} for a good manual solution.
#' @param parameter The diversification parameter to be mapped onto the tree. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param discrete Logical. Whether to discretize the distribution of reconstructed rates into bins
#' @param breaks A numeric vector of cut points for binning the rates. Passed internally to \code{cut}. The function checks whether 0 and max(rate) is in this vector and adds them if not. Scale appropriatelly if \code{plot_as_waiting_time=TRUE}.
#' @param colors Colors for the branches. If \code{discrete=TRUE} a vector of length matching the number of bins produced by the breaks. If \code{discrete=FALSE}, a vector of length two giving the min/max of a color gradient.
#' @param plot_as_waiting_time Whether to plot the rates (FALSE, default) or their inverse (waiting times)
#' @param tree_layout A layout for the tree. Available options are 'rectangular' (default), 'slanted', 'circular', 'fan' and 'radial'.
#' @param tree_direction 'right' (default), 'left', 'up', or 'down' for rectangular and slanted tree layouts
#' @param time_axis_ticks numeric giving the number of ticks for the time axis (default=10)
#' @param open_angle The degrees of empty space between the first and last tip. Only works for \code{tree_layout = 'fan'} and allows for a little more space around axis tick labels.
#'
#' @return A plot of the phylogeny with branches colored by hisse-inferred marginal ancestral states.
#'
#' @examples
#'
#'data("diatoms")
#'processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)
#'
#'h_rate_recon(
#'    processed_recon = processed_hisse,
#'    parameter = "extinction", discrete=FALSE)
#'
#'h_rate_recon(
#'    processed_recon = processed_hisse,
#'    parameter = "extinction", discrete=TRUE, breaks=c(0.3, 0.6, 1))
#'
#'@export

h_rate_recon <-
  function(processed_recon,
           show_tip_labels = FALSE,
           parameter = "turnover",
           discrete = FALSE,
           breaks = seq(0, 1, 0.2),
           colors = viridis(n=2),
           plot_as_waiting_time = FALSE,
           tree_layout = "rectangular",
           tree_direction = "right",
           time_axis_ticks = 10,
           open_angle = 10) {
    if (!tree_layout %in% c('rectangular', 'circular', 'slanted', 'fan', 'radial')) {
      stop("The selected tree layout is not supported.")
    }

    tree <- processed_recon$tree_data@phylo
    datas <- processed_recon$tree_data@data
    agemax <- tree %>% branching.times() %>% max()

    if (plot_as_waiting_time) {
      datas <- mutate(datas, "wanted" = 1 / !!as.name(parameter))
    } else {
      datas <- mutate(datas, "wanted" = !!as.name(parameter))
    }

    ggg <-
      ggtree(
        tr = tree,
        layout = tree_layout,
        size = 0.45,
        open.angle = 10
      ) +
      theme(
        # legend.direction = "horizontal",
        legend.position = "right",
        legend.key.size = unit(x = 0.5, units = "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank()
      )

    if (discrete) {
      max_rate <- datas %>% select(.data$wanted) %>% unlist %>% unname %>% max
      if (!0 %in% breaks) {
        breaks <- c(0, breaks)
      }
      if (all(max_rate > breaks)) {
        breaks <- c(breaks, max_rate)
      }
      message("Cutting distribution of rate with these breaks:\n")
      print(breaks)
      if (length(breaks != length(colors))) {
        message("The number of supplied colors does not match the number of breaks. Choosing different colors.")
        colors <- viridis(n = length(breaks))
      }
      param <- datas %>%
        select(.data$wanted) %>%
        unlist %>%
        unname %>%
        cut(breaks = breaks)
      ggg <-
        ggg + aes(color = param) +
        scale_color_manual(
          values=colors,
          name = parameter
        ) + guides(color = guide_legend(override.aes = list(size = 4)))
    } else {
      param <- datas %>% select(.data$wanted) %>% unlist %>% unname
      ggg <-
        ggg + aes(color = param) +
        scale_color_gradient(
          low = colors[1],
          high = colors[2],
          name = parameter
        )
    }

    ggg <-
      tree_flip(
        ggtree_object = ggg,
        show_tip_labels = show_tip_labels,
        tree_layout = tree_layout,
        tree_direction = tree_direction,
        time_axis_ticks = time_axis_ticks,
        agemax = agemax
      )

    return(ggg + theme(plot.margin = unit(rep(.1, 4), "in")))
  }

##### --- MuHiSSE functions ------------------------- #####

#' Extract transition rates from a multistate HiSSE model (MuSSE, MuHiSSE or CID)
#'
#'@description A helper function that formats \code{hisse::MuHiSSE}'s output into a matrix so transition rates between states can be easier to visualize.
#'
#'@param model_fit A muhisse fit object
#'@param hidden_states Logical whether the model contains hidden states
#'@param states A vector of state labels to replace hisse's 0/1 notation (see examples)
#'
#'@return The estimated transition rates formatted in a matrix
#'
#'@examples
#'
#'# The translation should follow the order of states in model_fit$solution.
#'# In this example, hisse's 00, 01, 10, and 11 states correspond to
#'# marine-plankton, marine-benthos, freshwater-plankton, and freshwater-benthos
#'
#'States <- c("mp", "mb", "fp", "fb")
#'
#'# muhisse model
#'data("diatoms")
#'m_transition_matrix(
#'  model_fit = diatoms$muhisse,
#'  hidden_states = TRUE,
#'  states = States)
#'
#'# CID8 model
#'# we have 8 hidden states, so the transition matrix is cumbersome
#'m_transition_matrix(
#'  model_fit = diatoms$cid8,
#'  hidden_states = TRUE,
#'  states = States)
#'
#'# musse model, no hidden states
#'m_transition_matrix(
#'  model_fit = diatoms$musse,
#'  hidden_states = FALSE,
#'  states = c("mp", "mb", "fp", "fb"))
#'
#'@export

m_transition_matrix <-
  function(model_fit,
           hidden_states = TRUE,
           states) {
    mat_size <- model_fit$trans.matrix %>% ncol
    all_states <-
      expand.grid(states, 1:(mat_size / 4)) %>%
      mutate(Name = paste(.data$Var1, .data$Var2, sep = "")) %>%
      select(.data$Name) %>% unlist %>% unname

    if (!hidden_states) {
      wanted <- c(
        "q00A_00A",
        "q01A_00A",
        "q10A_00A",
        "q11A_00A",
        "q00A_01A",
        "q01A_01A",
        "q10A_01A",
        "q11A_01A",
        "q00A_10A",
        "q01A_10A",
        "q10A_10A",
        "q11A_10A",
        "q00A_11A",
        "q01A_11A",
        "q10A_11A",
        "q11A_11A"
      )
    } else {
      tmp <-
        colnames(model_fit$trans.matrix) %>%
        str_replace(regex("\\("), "") %>%
        str_replace(regex("\\)"), "")
      wanted <- expand.grid(tmp, tmp) %>%
        mutate("Names" = paste("q", .data$Var1, "_", .data$Var2, sep = "")) %>%
        select(.data$Names) %>% unlist %>% unname
    }

    Rates <- numeric(length = length(wanted))
    names(Rates) <- wanted

    get_rate <- function(name, named_rates) {
      xx <- named_rates[which(names(named_rates) == name)]
      yy <- length(xx)
      if (yy > 0) {
        res <- xx
      } else {
        res <- NA
      }
      return(res)
    }

    # map(1:length(Rates), function(x) print(names(Rates)[x]))
    Rate_matrix <-
      map(1:length(Rates), function(x)
        get_rate(name = names(Rates)[x], named_rates = model_fit$solution)) %>%
      unlist %>% matrix(ncol = mat_size)

    Col_names <-
      colnames(model_fit$trans.matrix) %>%
      str_replace(regex("\\("), "") %>%
      str_replace(regex("\\)"), "")

    translation <- data.frame(Col_names, states_names = all_states)
    colnames(Rate_matrix) <-
      rownames(Rate_matrix) <- translation$states_names
    row_ind <- is.na(rep(states, (mat_size / 4)))
    col_ind <- is.na(rep(states, (mat_size / 4)))
    Rate_matrix <- Rate_matrix[!row_ind,!col_ind]
    if (!hidden_states) {
      colnames(Rate_matrix) <-
        rownames(Rate_matrix) <-
        str_remove(colnames(Rate_matrix), regex("\\d+"))
    }
    return(Rate_matrix)
  }

#' Extract estimated diversification rates from a multistate HiSSE model (MuSSE, MuHiSSE or CID)
#'
#'@description A helper function that formats \code{hisse::MuHiSSE}'s diversification parameter estimates into a data frame
#'
#'@param model_fit A muhisse fit object
#'@param states A vector of observed character states
#'
#'@return The estimated diversification rates formatted in a data.frame
#'
#'@examples
#'
#'States <- c("mp", "mb", "fp", "fb")
#'
#'# muhisse model
#'data("diatoms")
#'m_diversification_rates(model_fit = diatoms$muhisse, states=States)
#'
#'# CID8 model
#'m_diversification_rates(model_fit = diatoms$cid8, states=States)
#'
#'# musse model, no hidden states
#'m_diversification_rates(model_fit = diatoms$musse, states = States)
#'
#'@export

m_diversification_rates <- function(model_fit, states) {
  hidden_traits <-
    str_extract(string = colnames(model_fit$trans.matrix),
                pattern = regex("[A-Z]")) %>% unique
  hidden_traits <- as.numeric(as.factor(hidden_traits))
  hidden_traits <- ifelse(is.na(hidden_traits), 1, hidden_traits)
  num_hidden_traits <- length(unique(hidden_traits))

  num_states <- length(hidden_traits)
  num_states <- num_states * 4

  tur <-
    tibble("par_name" = names(model_fit$solution),
           "par_value" = model_fit$solution) %>%
    filter(stringr::str_detect(.data$par_name, "turnover")) %>%
    slice(1:num_states) %>%
    select(.data$par_value) %>% unlist %>% unname

  eps <-
    tibble("par_name" = names(model_fit$solution),
           "par_value" = model_fit$solution) %>%
    filter(stringr::str_detect(.data$par_name, "eps")) %>%
    slice(1:num_states) %>%
    select(.data$par_value) %>% unlist %>% unname

  States <- expand.grid(states, hidden_traits) %>%
    mutate("Name" = paste(.data$Var1, .data$Var2, sep = "")) %>%
    select(.data$Name) %>% unlist %>% unname

  res <-
    tibble(
      "State" = States,
      "Turnover" = tur,
      "Extinction_fraction" = eps
    )

  na_index <-
    rep(states, num_hidden_traits) %>% is.na() %>% as.numeric()
  res <-
    res %>% mutate(na_index = na_index) %>% filter(na_index == 0) %>% select(-na_index)
  res_len <- nrow(res)
  res <- res %>%
    mutate(Speciation =
             ParameterTransform(
               .data$Turnover,
               .data$Extinction_fraction)[1:res_len]
           ) %>%
    mutate(Extinction =
             ParameterTransform(
               .data$Turnover,
               .data$Extinction_fraction)[(res_len +1):(res_len + res_len)]
           ) %>%
    mutate(Net_diversification = .data$Speciation - .data$Extinction)

  return(res)
}

#' Extract transition and diversification rates from a multistate HiSSE model (MuSSE, MuHiSSE or CID)
#'
#'@description This is a function that wraps \code{m_transition_matrix} and \code{m_diversification_rates}
#' so rates can be extracted in one step
#'
#'@param model_fit A muhisse.fit object
#'@param hidden_traits Logical indicating whether the model has hidden states
#'@param character_states A translation for the character states in the model in the order 00, 01, 10, 11
#'
#'@return A list with two components, the transition rates formatted as a rate matrix and a data frame of diversification rates
#'
#'@examples
#'
#'States <- c("mp", "mb", "fp", "fb")
#'
#'# muhisse model
#'data("diatoms")
#'m_collect_rates(model_fit = diatoms$muhisse, hidden_traits=TRUE, character_states=States)
#'
#'@export

m_collect_rates <-
  function(model_fit,
           hidden_traits,
           character_states) {
    trans_rates <-
      m_transition_matrix(model_fit = model_fit,
                          hidden_states = hidden_traits,
                          states = character_states)
    div_rates <-
      m_diversification_rates(model_fit = model_fit, states = character_states)

    return(list(
      Transition_rates = trans_rates,
      Diversification_rates = div_rates
    ))
  }

#' Prepare MuHiSSE marginal reconstruction for plotting
#'
#'@description A function to process the output from \code{hisse::MarginReconMuHiSSE}. First we get model averaged rates for tips and nodes of the phylogeny and probabilities from the ancestral reconstruction of character states. Then convert the tree from \code{phylo} format into a \code{treedata} class and assign the data.frame of mode-averaged rates to the \code{@data} slot. We output the node and tip rates tables and the reformated tree object for downstream plotting with \code{ggtree} and \code{ggplot}.
#'
#'@param muhisse_recon An object produced by \code{hisse::MarginReconMuHiSSE}, or a list of such objects over which we can model-average using \code{hisse::GetModelAveRates}
#'
#'@details Internally, the probabilities for individual states (00, 01, 10, 11) are converted into \code{prob_0x} giving the probability for the first character being 0, and \code{prob_x0} giving the probability of the second character being 0. This way we can plot in terms of a 'focal character', i.e., we can plot the rate on the y axis and the probability of being 0 for the 'focal character' on the x axis while color-coding by the states of both characters.
#'
#'#'@return A list with three components: \itemize{
#'  \item \code{tree_data} the phylogeny with associated data object with reconstructed ancestral, model-averaged node and tip states and rates
#'  \item \code{tip_rates/node_rates} a tibble of reconstructed ancestral, model-averaged tip/node states and rates.
#'  }
#'
#'
#'@examples
#'
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'processed_muhisse$tip_rates
#'processed_muhisse$node_rates
#'processed_muhisse$tree_data
#'
#'@export

m_process_recon <- function(muhisse_recon) {
  tip.rates <-
    GetModelAveRates(x = muhisse_recon, type = "tips") %>%
    as_tibble() %>%
    mutate("prob_0x" = .data$state.00 + .data$state.01,
           "prob_x0" = .data$state.00 + .data$state.10)
  colnames(tip.rates)[1] <- "id"

  nod.rates <-
    GetModelAveRates(x = muhisse_recon, type = "nodes") %>%
    as_tibble() %>%
    mutate("prob_0x" = .data$state.00 + .data$state.01,
           "prob_x0" = .data$state.00 + .data$state.10)

  nod.rates$id <- as.character(nod.rates$id)
  tip.rates$id <- as.character(tip.rates$id)
  both.rates <- bind_rows(tip.rates, nod.rates)

  if (class(muhisse_recon) == "list") {
    tree <- as.treedata(muhisse_recon[[1]]$phy)
  }

  if (class(muhisse_recon) == "muhisse.states") {
    tree <- as.treedata(muhisse_recon$phy)
  }

  tree@data <- both.rates

  return(list(
    "tree_data" = tree,
    "tip_rates" = tip.rates,
    "node_rates" = nod.rates
  ))
}

#' Plot diversification rates estimated by a HiSSE model with means and standard deviations across tips and a two-dimensional colorplane for color
#'
#' @description  A function to plot a jittered scatterplot of (model-averaged) diversification rates in the alternative states.
#'
#' @param processed_recon An object produced with \code{m_process_recon}
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param focal_character Specifies the x axis. Either \code{prob_0x} to plot the probability of state 0 for the first character, or \code{prob_x0} to plot the probability for state 0 for the second character.
#' @param focal_character_label Label for the x axis of the scatterplot and two-dimensional color gradient. This should match the focal probability.
#' @param second_character_label Label for the y axis of the scatterplot and two-dimensional color gradient.
#' @param colors A vector of three colors in the order: (1) zero color (color when the two traits are in state 0), (2) horizontal_color (color to interpolate towards state 1 of the focal character) and (2) vertical_color (color to interpolate towards state 1 of the second character). See \code{?colorplaner::color_projections} for details.
#' @param plot_as_waiting_time Logical, whether to convert the rate to waiting time (1/rate)
#'
#' @return A scatterplot with focal probability (0 or 1) on the x axis and the chosen diversification parameter on the y axis with means and error bars (mean +/- SD) for each state color coded with in two-dimensional colorplane.
#'
#' @examples
#'
#'library("colorplaner")
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'m_scatterplot_cp(
#'  processed_recon = processed_muhisse,
#'  parameter = "turnover",
#'  focal_character = "prob_0x",
#'  focal_character_label = "p(mar)",
#'  second_character_label = "p(pla)",
#'  colors = c("#21908CFF", "#440154FF", "#FDE725FF"),
#'  plot_as_waiting_time = TRUE) +
#'  labs(y="Net turnover\n(waiting time in millions of years)")
#'
#'@export

m_scatterplot_cp <-
  function(processed_recon,
           parameter = "turnover",
           focal_character = c("prob_0x", "prob_x0"),
           focal_character_label,
           second_character_label,
           colors,
           plot_as_waiting_time = FALSE) {
    if (plot_as_waiting_time) {
      tip_rates <- processed_recon$tip_rates %>%
        mutate(wanted = 1 / !!as.name(parameter))
    } else
      tip_rates <- processed_recon$tip_rates %>%
        mutate(wanted = !!as.name(parameter))

    max_rate <- tip_rates %>% select(.data$wanted) %>% unlist %>% max

    sum_tip_rates <- tip_rates %>%
      mutate(both_prob = interaction(.data$prob_0x, .data$prob_x0)) %>%
      group_by(.data$both_prob) %>%
      select(.data$both_prob, .data$wanted) %>%
      summarise_at(.vars = vars("wanted"),
                   .funs = list("MN" = mean, "SD" = sd)) %>%
      mutate(LB = .data$MN - .data$SD, UB = .data$MN + .data$SD, wanted=.data$MN) %>%
      mutate(prob_0x = str_extract(.data$both_prob, regex("^\\d+")) %>% as.numeric()) %>%
      mutate(prob_x0 = str_extract(.data$both_prob, regex("\\d+$")) %>% as.numeric())

    if (focal_character == "prob_0x") {
      sum_tip_rates <-
        mutate(sum_tip_rates, focal_character = factor(.data$prob_0x))
      print(sum_tip_rates)
      sss <-
        ggplot(
          tip_rates,
          aes(
            x = factor(!!as.name(focal_character)),
            y = .data$wanted,
            color = .data$prob_0x,
            color2 = .data$prob_x0
          )
        )
      nudgex <- c(-0.3,-0.3, 0.3, 0.3)
    }

    if (focal_character == "prob_x0") {
      sum_tip_rates <-
        mutate(sum_tip_rates, focal_character = factor(.data$prob_x0))
      print(sum_tip_rates)
      sss <-
        ggplot(
          tip_rates,
          aes(
            x = factor(!!as.name(focal_character)),
            y = .data$wanted,
            color = .data$prob_x0,
            color2 = .data$prob_0x
          )
        )
      nudgex <- c(-0.3, 0.3,-0.3, 0.3)
    }

    if (nrow(sum_tip_rates) == 3) {
      nudgex <- rep(.3, 3)
    }

    sss <- sss +
      geom_point(
        alpha = .7,
        size = 2,
        position = position_jitter(width = .15, height = 0)
      ) +
      geom_errorbar(
        data = sum_tip_rates,
        aes(
          x = .data$focal_character,
          ymin = .data$LB,
          ymax = .data$UB
        ),
        position = position_nudge(x = nudgex, y = 0),
        inherit.aes = TRUE,
        width = .04
      ) +
      geom_point(
        data = sum_tip_rates,
        aes(x = .data$focal_character, y = .data$MN),
        position = position_nudge(x = nudgex, y = 0),
        inherit.aes = TRUE,
        pch = 21,
        stroke = 1,
        size = 2
      ) +
      scale_color_colorplane(
        color_projection = interpolate_projection,
        zero_color = colors[1],
        horizontal_color = colors[2],
        vertical_color = colors[3],
        axis_title = focal_character_label,
        axis_title_y = second_character_label,
        breaks = c(0, 0.5, 1),
        breaks_y = c(0, 0.5, 1),
        labels = as.character(c(0, 0.5, 1)),
        labels_y = as.character(c(0, 0.5, 1))
      ) +
      scale_y_continuous(breaks = pretty(x = c(0, max_rate), n = 10)) +
      theme_classic() +
      theme(legend.position = "right",
            legend.key.size = unit(x = .6, units = "cm")) +
      labs(x = focal_character_label, y = parameter) +
      theme_xy
    return(sss)
  }


#' Plot diversification rates estimated by a MuHiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a ridgeline of (model-averaged) diversification rates in the alternative states.
#'
#' @param processed_recon An object produced by \code{m_process_recon}
#' @param states_names Translation for the character states in the order 00, 01, 10, 11 (if wanted)
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param fill_colors Colors for the density polygons
#' @param line_colors Colors for the lines and points
#' @param plot_as_waiting_time Whether to plot the rates (FALSE, default) or their inverse (waiting times)
#'
#' @return A ridgeline plot of tip-associated rates (possibly model averaged), with carpet with the individual points and summary statistics (mean +/- sd)
#'
#' @examples
#'
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'
#'m_ridgelines(
#'  processed_recon = processed_muhisse,
#'  parameter = "extinction",
#'  line_colors = c("#000004FF", "#57157EFF", "#C43C75FF", "#FE9F6DFF"))
#'
#'@export

m_ridgelines <- function(processed_recon,
                         states_names = c("00", "01", "10", "11"),
                         parameter = "turnover",
                         plot_as_waiting_time = FALSE,
                         fill_colors = rep(NA, 4),
                         line_colors = viridis(n = 4)) {
  message(
    "Recoding and renaming character states. The elements 1:4 of the vector `character_states_names` are assumed to match the states 00, 01, 10, 11.\n"
  )

  ss <-
    m_prep_df(
      processed_recon = processed_recon,
      states_names = states_names,
      parameter = parameter
    )

  message("Summarising grouped by character state\n")
  wanted <- as.name(parameter)

  ss_var <-
    ss %>% group_by(.data$four_state) %>% select(.data$four_state, .data$wanted) %>% distinct()
  if (nrow(ss_var) == length(states_names)) {
    print(ss_var)
    stop(
      "Looks like there is no variation within the observed states. Probably because the model has no hidden states (e.g. MuSSE?). No point in calculating densities for point estimates."
    )
  }

  if (plot_as_waiting_time) {
    ss <- mutate(ss, wanted = 1 / !!wanted)
    ss.sum <- ss %>% group_by(.data$four_state) %>%
      summ()
  } else {
    ss <- mutate(ss, wanted = !!wanted)
    ss.sum <- ss %>% group_by(.data$four_state) %>%
      summ()
  }
  max_rate <-
    ss %>% select(.data$wanted) %>% top_n(1, wt = .data$wanted) %>% unlist %>% unname %>% unique
  # print(max_rate)
  print(ss.sum)

  message("\nPlotting\n\n")

  ggplot(
    data = ss,
    aes(
      x = .data$wanted,
      y = .data$four_state,
      fill = .data$four_state,
      colour = .data$four_state
    )
  ) +
    geom_density_ridges(
      alpha = 0.75,
      size = 0.75,
      jittered_points = TRUE,
      scale = .95,
      rel_min_height = .01,
      point_shape = "|",
      point_size = 1,
      position = position_nudge(y = rep(-.2, 4))
    ) +
    geom_errorbarh(
      data = ss.sum,
      position = position_nudge(y = rep(-0.3, 4)),
      aes(
        xmin = .data$Mean - .data$SD,
        xmax = .data$Mean + .data$SD,
        y = .data$four_state,
        colour = .data$four_state
      ),
      height = 0.05,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = ss.sum,
      pch = 21,
      position = position_nudge(y = rep(-0.3, 4)),
      aes(
        y = .data$four_state,
        x = .data$Mean,
        colour = .data$four_state
      ),
      size = 3,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(breaks = pretty(x = c(0, max_rate), n = 10)) +
    scale_fill_manual(name="", values = fill_colors) +
    scale_colour_manual(name="", values = line_colors) +
    labs(y = "", x = parameter) +
    theme_classic() +
    theme(legend.position = "none") +
    theme_xy
}

#' Plot diversification rates estimated by a MuHiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a jittered scatterplot of (model-averaged) diversification rates and summary statistics in the alternative states.
#'
#' @param processed_recon An object produced by \code{m_process_recon}
#' @param states_names Translation for the character states in the order 00, 01, 10, 11 (if wanted)
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param colors Colors for the lines and points
#' @param plot_as_waiting_time Whether to plot the rates (FALSE, default) or their inverse (waiting times)
#'
#' @return A jittered scatterplot of tip-associated rates (possibly model averaged) with summary statistics (mean +/- sd)
#'
#' @examples
#'
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'
#'m_scatterplot(
#'  processed_recon = processed_muhisse,
#'  parameter = "extinction",
#'  colors = c("#000004FF", "#57157EFF", "#C43C75FF", "#FE9F6DFF"))
#'
#'@export

m_scatterplot <-
  function(processed_recon,
           states_names = c("00", "01", "10", "11"),
           parameter = "turnover",
           colors = viridis(n = 4, option = 2, end = 0.7),
           plot_as_waiting_time = FALSE) {
    message(
      "Recoding and renaming character states. The elements 1:4 of the vector `character_states_names` are assumed to match the states 00, 01, 10, 11.\n"
    )

    ss <-
      m_prep_df(
        processed_recon = processed_recon,
        states_names = states_names,
        parameter = parameter
      )

    message("Summarising grouped by character state\n")
    wanted <- as.name(parameter)

    if (plot_as_waiting_time) {
      ss <- mutate(ss, wanted = 1 / !!wanted)
      ss.sum <- ss %>% group_by(.data$four_state) %>%
        summ()
    } else {
      ss <- mutate(ss, wanted = !!wanted)
      ss.sum <- ss %>% group_by(.data$four_state) %>%
        summ()
    }
    max_rate <-
      ss %>% select(.data$wanted) %>% top_n(1, wt = .data$wanted) %>% unlist %>% unname
    print(select(ss.sum, -.data$wanted))

    message("\nPlotting\n")

    pl <-
      ggplot(data = ss,
             aes(
               x = .data$four_state,
               y = .data$wanted,
               fill = .data$four_state
             )) +
      geom_point(
        pch = 21,
        position = position_jitter(width = 0.2, height = 0),
        size = 2,
        stroke = .3
      ) +
      geom_errorbar(
        inherit.aes = FALSE,
        data = ss.sum,
        aes(
          x = .data$four_state,
          colour = .data$four_state,
          ymin = .data$Mean - .data$SD,
          ymax = .data$Mean + .data$SD
        ),
        width = .05,
        position = position_nudge(x = .35)
      ) +
      geom_point(
        inherit.aes = FALSE,
        data = ss.sum,
        aes(
          x = .data$four_state,
          y = .data$Mean,
          colour = .data$four_state
        ),
        pch = 21,
        fill = "white",
        size = 3,
        stroke = 2,
        position = position_nudge(x = .35)
      ) +
      scale_fill_manual(name = "", values = colors) +
      scale_colour_manual(name = "", values = colors) +
      scale_y_continuous(breaks = pretty(x = c(0, max_rate), n = 10)) +
      theme_classic() +
      labs(x = "",
           y = parameter) +
      theme(legend.position = "top")+
      theme_xy
    return(pl)
  }

#' Plot diversification rates estimated by a MuHiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a stacked dotplot of (model-averaged) diversification rates and summary statistics in the alternative states.
#'
#' @param processed_recon An object produced by \code{m_process_recon}
#' @param states_names Translation for the character states in the order 00, 01, 10, 11 (if wanted)
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param colors Colors for the lines and points
#' @param plot_as_waiting_time Whether to plot the rates (FALSE, default) or their inverse (waiting times)
#' @param bin_width The width of bins for the dotplot. Treat this as any histogram. Testing several different bin width values is recommended.
#'
#' @return A stacked dotplot of tip-associated rates (possibly model averaged) with summary statistics (mean +/- sd)
#'
#' @examples
#'
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'
#'m_dotplot(
#'  processed_recon = processed_muhisse,
#'  parameter = "extinction",
#'  colors = c("#000004FF", "#57157EFF", "#C43C75FF", "#FE9F6DFF"))
#'
#'@export

m_dotplot <-
  function(processed_recon,
           states_names = c("00", "01", "10", "11"),
           parameter = "turnover",
           colors = viridis(n = 4, option = 2, end = 0.7),
           plot_as_waiting_time = FALSE,
           bin_width = 0.1) {
    message(
      "Recoding and renaming character states. The elements 1:4 of the vector `character_states_names` are assumed to match the states 00, 01, 10, 11.\n\n"
    )

    ss <-
      m_prep_df(
        processed_recon = processed_recon,
        states_names = states_names,
        parameter = parameter
      )

    message("Summarising grouped by character state\n\n")
    wanted <- as.name(parameter)

    ss_var <-
      ss %>% group_by(.data$four_state) %>% select(.data$four_state, .data$wanted) %>% distinct()
    if (nrow(ss_var) == length(states_names)) {
      print(ss_var)
      stop(
        "Looks like there is no variation within the observed states. Probably because the model has no hidden states (e.g. MuSSE?). No point in calculating histograms for point estimates."
      )
    }

    if (plot_as_waiting_time) {
      ss <- mutate(ss, wanted = 1 / !!wanted)
      ss.sum <- ss %>% group_by(.data$four_state) %>%
        summ()
    } else {
      ss <- mutate(ss, wanted = !!wanted)
      ss.sum <- ss %>% group_by(.data$four_state) %>%
        summ()
    }
    max_rate <-
      ss %>% select(.data$wanted) %>% top_n(1, wt = .data$wanted) %>% unlist %>% unname
    # print(ss)
    print(select(ss.sum,-.data$wanted))

    message("Plotting")

    pl <- ggplot(data = ss,
                 aes(
                   x = .data$four_state,
                   y = .data$wanted,
                   fill = .data$four_state
                 )) +
      geom_dotplot(
        alpha = .75,
        colour = "white",
        binwidth = bin_width,
        dotsize = 1.2,
        stackratio = 0.5,
        binpositions = "all",
        binaxis = "y",
        stackdir = "up"
      ) +
      geom_errorbar(
        inherit.aes = FALSE,
        data = ss.sum,
        aes(
          x = .data$four_state,
          colour = .data$four_state,
          ymin = .data$Mean - .data$SD,
          ymax = .data$Mean + .data$SD
        ),
        width = .05,
        position = position_nudge(x = -.1)
      ) +
      geom_point(
        inherit.aes = FALSE,
        data = ss.sum,
        aes(
          x = .data$four_state,
          y = .data$Mean,
          colour = .data$four_state
        ),
        pch = 21,
        fill = "white",
        size = 3,
        stroke = 2,
        position = position_nudge(x = -.1)
      ) +
      scale_fill_manual(name = "", values = colors) +
      scale_colour_manual(name = "", values = colors) +
      scale_y_continuous(breaks = pretty(x = c(0, max_rate), n = 10)) +
      labs(x = "", y = parameter) +
      theme_classic() +
      theme(legend.position = "none") +
      theme_xy
    return(pl)
  }

#' Plot MuHiSSE model-averaged marginal ancestral state reconstruction for the trait with a two-dimensional colorplane
#'
#' @description A function to plot a MuHiSSE (model-averaged) marginal ancestral reconstruction for the trait data.
#'
#' @param processed_recon An object produced by \code{m_process_recon}
#' @param show_tip_labels Logical, whether to plot tip labels. Default is FALSE because it is difficult to plot legible tip labels for larger trees common in this type of analysis. See \code{?m_trait_recon} for a good manual solution.
#' @param tree_layout A layout for the tree. Available options are 'rectangular' (default), 'slanted', 'circular', 'fan' and 'radial'.
#' @param tree_direction 'right' (default), 'left', 'up', or 'down' for rectangular and slanted tree layouts
#' @param time_axis_ticks numeric giving the number of ticks for the time axis (default=10). Passed on to \code{pretty} internally, so the number of plotted ticks might not be exactly the same.
#' @param open_angle The degrees of empty space between the first and last tip. Only works for \code{tree_layout = 'fan'} and allows for a little more space around axis tick labels.
#' @param focal_character Specifies the x axis of the two-dimensional color key. Either \code{prob_0x} to plot the probability of state 0 for the first character, or \code{prob_x0} to plot the probability for state 0 for the second character.
#' @param focal_character_label Label for the x axis of the two-dimensional color gradient. This should match the focal probability.
#' @param second_character_label Label for the y axis of the two-dimensional color gradient.
#' @param colors A vector of three colors in the order: (1) zero color (color when the two traits are in state 0), (2) horizontal_color (color to interpolate towards state 1 of the focal character) and (2) vertical_color (color to interpolate towards state 1 of the second character). See \code{?colorplaner::color_projections} for details.
#'
#' @return A plot of the phylogeny with branches colored by muhisse-inferred marginal ancestral states.
#'
#' @examples
#'
#'library("colorplaner")
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'
#'m_trait_recon_cp(
#'  processed_recon = processed_muhisse,
#'  tree_layout = "fan",
#'  focal_character = "prob_0x",
#'  focal_character_label = "p(marine)",
#'  second_character_label = "p(plankton)",
#'  colors = c("#21908CFF", "#440154FF", "#FDE725FF")
#'  )
#'
#'@export

m_trait_recon_cp <-
  function(processed_recon,
           show_tip_labels = FALSE,
           tree_layout = "rectangular",
           tree_direction = "right",
           time_axis_ticks = 10,
           open_angle = 10,
           focal_character = c("prob_0x", "prob_x0"),
           focal_character_label,
           second_character_label,
           colors = viridis(9)[c(5, 1, 9)]) {
    if (!tree_layout %in% c('rectangular', 'circular', 'slanted', 'fan', 'radial')) {
      stop("The selected tree layout is not supported.")
    }

    tree <- processed_recon$tree_data@phylo
    datas <- processed_recon$tree_data@data
    agemax <- tree %>% branching.times() %>% max()

    ggg <-
      ggtree(
        tr = tree,
        layout = tree_layout,
        size = .45,
        open.angle = open_angle
      ) +
      theme(
        legend.direction = "horizontal",
        legend.position = "right",
        legend.key.size = unit(x = .5, units = "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(),
        plot.background = element_blank()
      ) +
      scale_color_colorplane(
        color_projection = interpolate_projection,
        zero_color = colors[1],
        horizontal_color = colors[2],
        vertical_color = colors[3],
        axis_title = focal_character_label,
        axis_title_y = second_character_label,
        breaks = c(0, 0.5, 1),
        breaks_y = c(0, 0.5, 1),
        labels = as.character(c(0, 0.5, 1)),
        labels_y = as.character(c(0, 0.5, 1))
      )

    if (focal_character == "prob_0x") {
      ggg <- ggg + aes(color = datas$prob_0x, color2 = datas$prob_x0)
    }

    if (focal_character == "prob_x0") {
      ggg <- ggg + aes(color = datas$prob_x0, color2 = datas$prob_0x)
    }

    ggg <-
      tree_flip(
        ggtree_object = ggg,
        show_tip_labels = show_tip_labels,
        tree_layout = tree_layout,
        tree_direction = tree_direction,
        time_axis_ticks = time_axis_ticks,
        agemax = agemax
      )

    return(ggg + theme(plot.margin = unit(rep(.1, 4), "in")))
  }


#' Plot MuHiSSE model-averaged marginal ancestral state reconstruction for the trait with discretized probabilities
#'
#' @description A function to plot a MuHiSSE (model-averaged) marginal ancestral reconstruction for the trait data.
#'
#' @param processed_recon An object produced by \code{m_process_recon}
#' @param show_tip_labels Logical, whether to plot tip labels. Default is FALSE because it is difficult to plot legible tip labels for larger trees common in this type of analysis. See \code{?m_trait_recon} for a good manual solution.
#' @param tree_layout A layout for the tree. Available options are 'rectangular' (default), 'slanted', 'circular', 'fan' and 'radial'.
#' @param tree_direction 'right' (default), 'left', 'up', or 'down' for rectangular and slanted tree layouts
#' @param time_axis_ticks numeric giving the number of ticks for the time axis (default=10). Passed on to \code{pretty} internally, so the number of plotted ticks might not be exactly the same.
#' @param open_angle The degrees of empty space between the first and last tip. Only works for \code{tree_layout = 'fan'} and allows for a little more space around axis tick labels.
#' @param colors A vector of colors for the character states. Note that the number of categories depends on the cutoff value. Internally this vector of colors is subset to \code{colors[1:number of states]}, so some tinkering might be necessary to get colors that work well together.
#' @param cutoff A vector of length 2 giving the values for discretizing the reconstructed probabilities for the first and second character. Categories "0", "uncertain", or "1" are assigned to probabilities \code{>= 1-cutoff}, \code{1-cutoff:cutoff}, \code{<=cutoff}. A cutoff=0.5 will result with four categories, essentially ignoring the uncertainty in ancestral state reconstruction. A cutoff < 0.5 will produce a category "uncertain" for probabilities in the range \code{1-cutoff:cutoff}. The function will print a count of the categories after discretizing, so it is possible to adjust the cutoff values to avoid plotting of categories with too few nodes
#' @param states_of_first_character Translation for the character states
#' @param states_of_second_character Translation for the character states
#'
#' @return A plot of the phylogeny with branches colored by muhisse-inferred marginal ancestral states.
#'
#' @examples
#'
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'cols = c("#00204DFF", "#1B3B6DFF", "#4E576CFF", "#727274FF", "#958F78FF", "#BCAF6FFF", "#E7D159FF")
#'# eight categories after binning with a cutoff of 0.2
#'m_trait_recon(
#'  processed_recon = processed_muhisse,
#'  cutoff = c(.2, .2),
#'  states_of_first_character = c("marine", "freshwater"),
#'  states_of_second_character = c("plankton", "benthos"),
#'  tree_layout = "radial",
#'  colors = cols)
#'
#'# three of these eight have < 3 nodes, so we could try to avoid plotting some of them
#'# adjust the cutoff for the second variable
#'
#'m_trait_recon(
#'  processed_recon = processed_muhisse,
#'  cutoff = c(.2, .3),
#'  states_of_first_character = c("marine", "freshwater"),
#'  states_of_second_character = c("plankton", "benthos"),
#'  tree_layout = "radial",
#'  colors = cols)
#'
#'# ignoring uncertainty
#'m_trait_recon(
#'  processed_recon = processed_muhisse,
#'  cutoff = c(.5, .5),
#'  states_of_first_character = c("marine", "freshwater"),
#'  states_of_second_character = c("plankton", "benthos"),
#'  tree_layout = "radial",
#'  colors = cols)
#'
#'
#'
#'#'# Plotting tip labels
#'# for better control over font size, justification, alignment, and offset
#'# use show_tip_labels= FALSE and add the tip labels manually
#'#
#'# for example
#'
#'require("ape")
#'data("diatoms")
#'x <- m_process_recon(diatoms$`3state_musse_recon`)
#'
#'# the defaults work poorly
#'m_trait_recon(processed_recon = x,
#'  show_tip_labels= TRUE,
#'  states_of_first_character = c("0", "1"),
#'  states_of_second_character = c("0", "1"))
#'
#'# manual setup is better
#'basic_plot <-
#'  m_trait_recon(
#'  processed_recon = x,
#'  show_tip_labels= FALSE,
#'  states_of_first_character = c("0", "1"),
#'  states_of_second_character = c("0", "1"),
#'  time_axis_ticks = 8
#')
#'
#'# extract the data component from the ggplot object
#'# we'll need it to order the tip labels
#'plot_data <- basic_plot$data %>%
#'  filter(isTip == TRUE) %>% # we need the tips only, no nodes
#'  arrange(y) # sort numerically
#'
#'basic_plot +
#'  #now add an axis to the right
#'  scale_y_continuous(position = "right",
#'                     breaks = plot_data$y,
#'                     labels = plot_data$label) +
#'  # enable the axis text, but disable the line and ticks marks
#'  # move the legend out of the way
#'  theme(
#'    legend.position = "top",
#'    axis.ticks.y = element_blank(),
#'    axis.line.y = element_blank(),
#'    axis.text.y = element_text(face = "italic")
#'  )
#'
#'@export

m_trait_recon <-
  function(processed_recon,
           show_tip_labels = FALSE,
           cutoff = c(.2, .2),
           states_of_first_character,
           states_of_second_character,
           tree_layout = "rectangular",
           tree_direction = "right",
           time_axis_ticks = 10,
           open_angle = 10,
           colors = viridis(n = 9)) {
    if (!tree_layout %in% c('rectangular', 'circular', 'slanted', 'fan', 'radial')) {
      stop("The selected tree layout is not supported.")
    }

    tree <- processed_recon$tree_data@phylo
    agemax <- tree %>% branching.times() %>% max()

    ss <- processed_recon$tree_data@data %>%
      mutate(
        prob_0x_named =
          case_when(
            prob_0x >= 1 - cutoff[1] ~ states_of_first_character[1],
            prob_0x <= cutoff[1] ~ states_of_first_character[2],
            TRUE ~ paste(
              states_of_first_character[1],
              "/",
              states_of_first_character[2],
              " uncertain",
              sep = ""
            )
          )
      ) %>%
      mutate(
        prob_x0_named =
          case_when(
            prob_x0 >= 1 - cutoff[2] ~ states_of_second_character[1],
            prob_x0 <= cutoff[2] ~ states_of_second_character[2],
            TRUE ~ paste(
              states_of_second_character[1],
              "/",
              states_of_second_character[2],
              " uncertain",
              sep = ""
            )
          )
      ) %>%
      mutate(four_state = paste(.data$prob_0x_named, .data$prob_x0_named, sep = "-"))

    message("Categories after discretizing with the provided cutoff:\n")
    ss.cnt <- ss %>%
      ungroup %>%
      group_by(.data$four_state) %>%
      add_tally() %>%
      select(.data$four_state, n) %>%
      distinct

    print(ss.cnt)

    nstat <- ss %>% select(.data$four_state) %>% distinct %>% nrow

    ggg <-
      ggtree(
        tr = tree,
        layout = tree_layout,
        size = .45,
        open.angle = open_angle,
        aes(color = ss$four_state)
      ) +
      scale_color_manual(name = "", values = colors[1:nstat]) +
      theme(
        legend.position = "right",
        legend.key.size = unit(x = .5, units = "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(),
        plot.background = element_blank()
      ) +
      guides(color = guide_legend(override.aes = list(size = 4)))

    ggg <-
      tree_flip(
        ggtree_object = ggg,
        show_tip_labels = show_tip_labels,
        tree_layout = tree_layout,
        tree_direction = tree_direction,
        time_axis_ticks = time_axis_ticks,
        agemax = agemax
      )

    return(ggg + theme(plot.margin = unit(rep(.1, 4), "in")))
  }

#' Plot MuHiSSE model-averaged marginal ancestral state reconstruction for diversification rates
#'
#' @description A function to plot a (model-averaged) marginal ancestral reconstruction for the estimated diversification rates.
#'
#' @param processed_recon An object produced by \code{h_process_recon}
#' @param show_tip_labels Logical, whether to plot tip labels. Default is FALSE because it is difficult to plot legible tip labels for larger trees common in this type of analysis. See \code{?m_trait_recon} for a good manual solution.
#' @param parameter The diversification parameter to be mapped onto the tree. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param discrete Logical. Whether to discretize the distribution of reconstructed rates into bins
#' @param breaks A numeric vector of cut points for binning the rates. Passed internally to \code{cut}. The function checks whether 0 and max(rate) is in this vector and adds them if not. Scale appropriatelly if \code{plot_as_waiting_time=TRUE}.
#' @param colors Colors for the branches. If \code{discrete=TRUE} a vector of length matching the number of bins produced by the breaks. If \code{discrete=FALSE}, a vector of length two giving the min/max of a color gradient.
#' @param plot_as_waiting_time Whether to plot the rates (FALSE, default) or their inverse (waiting times)
#' @param tree_layout A layout for the tree. Available options are 'rectangular' (default), 'slanted', 'circular', 'fan' and 'radial'.
#' @param tree_direction 'right' (default), 'left', 'up', or 'down' for rectangular and slanted tree layouts
#' @param time_axis_ticks numeric giving the number of ticks for the time axis (default=10)
#' @param open_angle The degrees of empty space between the first and last tip. Only works for \code{tree_layout = 'fan'} and allows for a little more space around axis tick labels.
#'
#' @return A plot of the phylogeny with branches colored by hisse-inferred marginal ancestral states.
#'
#' @examples
#'
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'
#'m_rate_recon(
#'    processed_recon = processed_muhisse,
#'    parameter = "extinction", discrete=FALSE)
#'
#'m_rate_recon(
#'    processed_recon = processed_muhisse,
#'    parameter = "extinction", discrete=TRUE, breaks=c(0.3, 0.6, 1))
#'
#'@export

m_rate_recon <-
  function(processed_recon,
           show_tip_labels = FALSE,
           parameter = "turnover",
           discrete = FALSE,
           breaks = seq(0.3, 0.6, 1),
           colors = c("red", "blue", "orange", "green"),
           plot_as_waiting_time = FALSE,
           tree_layout = "rectangular",
           tree_direction = "right",
           time_axis_ticks = 10,
           open_angle = 10) {
    if (!tree_layout %in% c('rectangular', 'circular', 'slanted', 'fan', 'radial')) {
      stop("The selected tree layout is not supported.")
    }

    tree <- processed_recon$tree_data@phylo
    datas <- processed_recon$tree_data@data
    agemax <- tree %>% branching.times() %>% max()

    if (plot_as_waiting_time) {
      datas <- mutate(datas, "wanted" = 1 / !!as.name(parameter))
    } else {
      datas <- mutate(datas, "wanted" = !!as.name(parameter))
    }

    ggg <-
      ggtree(
        tr = tree,
        layout = tree_layout,
        size = .45,
        open.angle = 10
      ) +
      theme(
        legend.position = "right",
        legend.key.size = unit(x = .5, units = "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank()
      )

    if (discrete) {
      max_rate <- datas %>% select(.data$wanted) %>% unlist %>% unname %>% max
      if (!0 %in% breaks) {
        breaks <- c(0, breaks)
      }
      if (all(max_rate > breaks)) {
        breaks <- c(breaks, max_rate)
      }
      message("Cutting distribution of rate with these breaks:\n")
      print(breaks)
      param <-
        datas %>%
        select(.data$wanted) %>%
        unlist %>%
        unname %>%
        cut(breaks = breaks)
      ggg <-
        ggg + aes(color = param) + scale_color_manual(
          values = colors,
          name = parameter
        ) + guides(color = guide_legend(override.aes = list(size = 4)))
    } else {
      param <- datas %>% select(.data$wanted) %>% unlist %>% unname
      ggg <-
        ggg + aes(color = param) + scale_color_gradient(
          low = colors[1],
          high = colors[2],
          name = parameter
        )
    }
    ggg <-
      tree_flip(
        ggtree_object = ggg,
        show_tip_labels = show_tip_labels,
        tree_layout = tree_layout,
        tree_direction = tree_direction,
        time_axis_ticks = time_axis_ticks,
        agemax = agemax
      )

    return(ggg + theme(plot.margin = unit(rep(.1, 4), "in")))
  }


##### --- Utility functions ------------------------- #####


tree_flip <- function(ggtree_object,
                      show_tip_labels,
                      tree_layout,
                      tree_direction,
                      time_axis_ticks,
                      agemax) {
  if (tree_layout %in% c("rectangular", "slanted")) {
    if (tree_direction == "up") {
      ggtree_object <- ggtree_object +
        coord_flip() +
        scale_x_continuous(
          expand = c(add = c(.0, .01)),
          breaks = pretty(c(0, agemax), n = time_axis_ticks),
          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))
        ) +
        theme(
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.text.y = element_text(size=12)
        ) +
        labs(x = "Time (Ma)")

      if (show_tip_labels) {
        plot_data <- ggtree_object$data %>%
          filter(.data$isTip == TRUE)
        ggtree_object <-
          ggtree_object +
          scale_y_continuous(
            expand = c(add = c(.02, .02)),
            position = "right",
            breaks = plot_data$y,
            labels = plot_data$label
          ) +
          theme(
            legend.position = "bottom",
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_text(
              angle = 90,
              hjust = 0,
              size = 3
            )
          )
      }
    }

    if (tree_direction == "down") {
      ggtree_object <- ggtree_object +
        coord_flip() +
        scale_x_continuous(
          trans = "reverse",
          expand = c(add = c(.0, .01)),
          breaks = pretty(c(0, agemax), n = time_axis_ticks),
          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))
        ) +
        theme(
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.text.y = element_text(size=12)
        ) +
        labs(x = "Time (Ma)")

      if (show_tip_labels) {
        plot_data <- ggtree_object$data %>%
          filter(.data$isTip == TRUE)
        ggtree_object <-
          ggtree_object +
          scale_y_continuous(
            expand = c(add = c(.02, .02)),
            position = "left",
            breaks = plot_data$y,
            labels = plot_data$label
          ) +
          theme(
            legend.position = "top",
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_text(
              angle = 270,
              hjust = 0,
              size = 3
            )
          )
      }
    }

    if (tree_direction == "left") {
      ggtree_object <- ggtree_object +
        scale_x_continuous(
          trans = "reverse",
          expand = c(add = c(.0, .01)),
          breaks = pretty(c(0, agemax), n = time_axis_ticks),
          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))
        ) +
        theme(
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(size=12)
        ) +
        labs(x = "Time (Ma)")

      if (show_tip_labels) {
        plot_data <- ggtree_object$data %>%
          filter(.data$isTip == TRUE)
        ggtree_object <-
          ggtree_object +
          scale_y_continuous(
            expand = c(add = c(.02, .02)),
            position = "left",
            breaks = plot_data$y,
            labels = plot_data$label
          ) +
          theme(
            legend.position = "top",
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.y = element_text(hjust = 1, size = 3)
          )
      }
    }

    if (tree_direction == "right") {
      ggtree_object <- ggtree_object +
        scale_x_continuous(
          expand = c(add = c(.0, .01)),
          breaks = pretty(c(0, agemax), n = time_axis_ticks),
          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))
        ) +
        theme(
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(size=12)
        ) +
        labs(x = "Time (Ma)")

      if (show_tip_labels) {
        plot_data <- ggtree_object$data %>%
          filter(.data$isTip == TRUE)
        ggtree_object <-
          ggtree_object +
          scale_y_continuous(
            expand = c(add = c(.02, .02)),
            position = "right",
            breaks = plot_data$y,
            labels = plot_data$label
          ) +
          theme(
            legend.position = "top",
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.y = element_text(hjust = 0, size = 3)
          )
      }
    }

  }

  if (tree_layout %in% c("circular", "fan", "radial")) {
    maxx <- ggtree_object$data %>%
      top_n(n = 1, wt = .data$x) %>%
      select(.data$x) %>%
      unlist %>%
      unname %>%
      unique
    maxx <- round(maxx, 1)
    ntip <-
      ggtree_object$data %>% filter(.data$isTip == TRUE) %>% nrow()
    ntip <- ntip + 10
    pretty_points <-
      maxx - c(maxx, pretty(c(maxx:0), n = time_axis_ticks))

    pp <- tibble(x = rev(pretty_points), y = 0) %>%
      filter(.data$x <= maxx, .data$x > 0) %>%
      mutate(label = rev(.data$x) - min(.data$x))

    ggtree_object <- ggtree_object +
      geom_vline(
        data = pp,
        aes(xintercept = .data$x),
        size = .2,
        color = "darkgrey"
      ) +
      geom_text(
        data = pp,
        aes(
          x = .data$x + 0.1,
          y = ntip + 2,
          label = .data$label
        ),
        size = 4,
        inherit.aes = FALSE
      )
    if (show_tip_labels) {
      ggtree_object <-
        ggtree_object + geom_tiplab2(inherit.aes = FALSE, size = .4)
    }
  }

  return(ggtree_object)
}

m_prep_df <-
  function(processed_recon,
           states_names,
           parameter) {
    ss <- processed_recon$tip_rates %>%
      mutate(
        what_state = paste(
          .data$state.00,
          .data$state.01,
          .data$state.10,
          .data$state.11,
          sep = ""
        )
      ) %>%
      mutate(four_state = factor(
        case_when(
          what_state == "1000" ~ states_names[1],
          what_state == "0100" ~ states_names[2],
          what_state == "0010" ~ states_names[3],
          what_state == "0001" ~ states_names[4]
        ),
        levels = states_names
      )) %>%
      mutate(wanted = !!as.name(parameter))

    return(ss)
  }


summ <- function(x) {
  y <- x %>% summarise_at(.vars = vars(.data$wanted),
                          .funs = list(
                            "Mean" = mean,
                            "Median" = median,
                            "SD" = sd
                          )) %>%
    mutate(wanted = .data$Mean) # to bypass warning in geom_errorbar
  return(y)
}

# not exported from hisse, not exported here
ParameterTransform <- function (x, y)
{
  speciation <- x / (1 + y)
  extinction <- (x * y) / (1 + y)
  return(c(speciation, extinction))
}


theme_xy <- theme(axis.text = element_text(size=12), axis.title = element_text(size=14))
