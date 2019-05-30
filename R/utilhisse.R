#' @name utilhisse
#' @docType package
#' @keywords package
#'
#' @import hisse dplyr ggplot2
#'
#' @importFrom stats median quantile sd
#' @importFrom stringr str_replace str_extract regex
#' @importFrom ggtree ggtree
#' @importFrom treeio as.treedata
#' @importFrom ggridges geom_density_ridges
#' @importFrom colorplaner scale_color_colorplane interpolate_projection
#' @importFrom viridisLite viridis
#' @importFrom viridis viridis
#' @importFrom ape branching.times
#' @importFrom purrr map
#'
NULL


##### --- TO DO ----------------------------------- #####

# 5. Support region functions. Making tables, plots, and contour plots?
# 6. g_ functions for HiGeoSSE
# 7. diagram functions
# 8. shiny app
# 9. test 3 state models

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
#' @param processed_hisse_recon An object produced by \code{h_process_recon}
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param x_label Label for the x axis. This is the binary trait whose assosiation with diversification is being tested.
#' @param states_names A character vector of length two giving the translation for states 0 and 1.
#' @param plot_as_waiting_time Logical. Whether to plot the rates or their inverse (waiting times)
#'
#' @return A jittered scatterplot of (model averaged) tip-associated rates.
#'
#' @examples
#'
#'data("diatoms")
#'processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)
#'hisse_rates_plot <- h_scatterplot(
#'  processed_hisse_recon=processed_hisse,
#'  parameter="turnover",
#'  x_label="habitat")
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
  function(processed_hisse_recon,
           parameter = "turnover",
           x_label = "",
           states_names = c("0", "1"),
           plot_as_waiting_time = FALSE) {
    tip.rates <- processed_hisse_recon$tip_rates
    tip.rates$f_state <- as.factor(tip.rates$state)

    if (plot_as_waiting_time) {
      tip.rates <- mutate(tip.rates, wanted = 1 / !!as.name(parameter))
    } else {
      tip.rates <- mutate(tip.rates, wanted = !!as.name(parameter))
    }

    tip.rates.sum <- tip.rates %>%
      group_by(f_state) %>%
      select(f_state, wanted) %>%
      summarise_if(
        .predicate = is.numeric,
        .funs = list("Mean" = mean,
                     "SD" = sd,
                     "Max" = max)) %>%
      mutate("wanted" = Mean) #to bypass warning in geom_errorbar

    result <-
      ggplot(data = tip.rates,
             aes(x = f_state,
                 y = wanted,
                 color = f_state)) +
      geom_point(alpha = .7,
                 size = 0.75,
                 position = position_jitter(width = .15)) +
      geom_errorbar(
        data = tip.rates.sum,
        width = .05,
        aes(
          x = f_state,
          ymax = Mean + SD,
          ymin = Mean - SD
        ),
        position = position_nudge(x = 0.3, y = 0)
      ) +
      geom_point(
        data = tip.rates.sum,
        size = 2.5,
        position = position_nudge(x = 0.3, y = 0),
        aes(x = f_state, y = Mean),
        pch = 21
      ) +
      scale_color_viridis_d(name = x_label,
                          end = 0.6) +
      scale_y_continuous(breaks = pretty(c(0, max(
        tip.rates.sum$Max
      )), n = 8)) +
      scale_x_discrete(breaks = c(0, 1), labels = states_names) +
      theme_classic() +
      theme(legend.position = "right",
            legend.key.size = unit(x = .6, units = "cm")) +
      labs(x = x_label, y = parameter)
    return(result)
  }

#' Plot diversification rates estimated by a HiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a dotplot of (model-averaged) diversification rates in the alternative states.
#'
#' @param processed_hisse_recon An object produced by \code{h_process_recon}
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param x_label Label for the x axis. This is the binary trait whose assosiation with diversification is being tested.
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
#'paint_cols <- c("orange", "violet")
#'
#' h_dotplot(
#'   processed_hisse_recon = processed_hisse,
#'   parameter = "turnover",
#'   states_names = c("Plankton", "Benthos"),
#'   bin_width = 0.2,
#'   colors=paint_cols[1:2],
#'   plot_as_waiting_time = TRUE
#' ) + labs(x = "", y = "waiting time (My)", title = "Turnover")
#'
#'# see ?h_plot_rates_states for examples for modifying the graph using ggplot2
#'
#'@export

h_dotplot <-
  function(processed_hisse_recon,
           parameter = "turnover",
           states_names = c("Marine", "Freshwater"),
           x_label = "",
           bin_width = 0.1,
           colors,
           plot_as_waiting_time = TRUE) {
    tip.rates <- processed_hisse_recon$tip_rates
    tip.rates$f_state <- as.factor(tip.rates$state)

    if (plot_as_waiting_time) {
      tip.rates <- mutate(tip.rates, "wanted" = 1 / !!as.name(parameter))
    } else {
      tip.rates <- mutate(tip.rates, "wanted" = !!as.name(parameter))
    }

    tip.rates.sum <- tip.rates %>%
      group_by(f_state) %>%
      select(f_state, wanted) %>%
      summarise_if(
        .predicate = is.numeric,
        .funs = list("Mean" = mean,
                     "SD" = sd,
                     "Max" =max)) %>%
      mutate("wanted" = Mean) #to bypass warning in geom_errorbar

    sss <-
      ggplot(data = tip.rates,
             aes(x = f_state,
                 y = wanted,
                 fill = f_state)) +
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
          x = f_state,
          ymax = Mean + SD,
          ymin = Mean - SD
        ),
        position = position_nudge(x = -0.1, y = 0)
      ) +
      geom_point(
        data = tip.rates.sum,
        size = 2.5,
        position = position_nudge(x = -0.1, y = 0),
        aes(x = f_state, y = Mean),
        pch = 21
      ) +
      scale_color_manual(values = colors, name = "") +
      scale_fill_manual(values = colors, name = "") +
      scale_x_discrete(breaks = c(0, 1), labels = states_names) +
      scale_y_continuous(breaks = pretty(c(0, max(
        tip.rates.sum$Max
      )), n = 8)) +
      theme_classic() +
      theme(
        legend.position = "none",
        legend.background = element_blank(),
        legend.key.size = unit(x = .6, units = "cm")
      ) +
      labs(x = x_label, y = parameter)
    return(sss)
  }

#' Plot diversification rates estimated by a HiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a ridgeline of (model-averaged) diversification rates in the alternative states. We can change the rate plotted on the y axis and modify the label for the x-axis (your binary 0/1 trait). Further modifications are straightforward with ggplot (see examples).
#'
#' @param processed_hisse_recon An object produced by \code{h_process_recon}
#' @param parameter The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param states_names The names for character states
#' @param colors Colors for the points in the two alternate states
#' @param plot_as_waiting_time Whether to plot the rates (FALSE, default) or their inverse (waiting times)
#'
#' @return A ridgeline plot of tip-associated rates (possibly model averaged).
#'
#' @examples
#'
#'data("diatoms")
#'processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)
#'paint_cols <- c("orange", "violet")
#'
#'h_ridgelines(
#'  processed_hisse_recon = processed_hisse,
#'  states_names = c("Plankton", "Benthos"),
#'  parameter = "extinction",
#'  colors = c("yellow", "red"))
#'
#'@export

h_ridgelines <- function(processed_hisse_recon,
                         parameter = "turnover",
                         states_names = c("Marine", "Freshwater"),
                         colors = c("yellow", "red"),
                         plot_as_waiting_time = FALSE) {
  tip.rates <- processed_hisse_recon$tip_rates
  tip.rates$f_state <-
    factor(ifelse(tip.rates$state == 0,
                  states_names[1],
                  states_names[2]),
           levels = states_names)

  if (plot_as_waiting_time) {
    tip.rates <- mutate(tip.rates, "wanted" = 1/!!as.name(parameter))
  } else {
    tip.rates <- mutate(tip.rates, "wanted" = !!as.name(parameter))
  }

  tip.rates.sum <- tip.rates %>%
    group_by(f_state) %>%
    select(f_state, wanted) %>%
    summarise_at(
      .vars = vars(wanted),
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
         aes(x = wanted,
             y = f_state,
             fill = f_state)) +
    geom_density_ridges(alpha = 0.75) +
    geom_point(
      data = tip.rates.sum,
      pch = 21,
      aes(y = c(1.2, 2.2), x = Mean),
      size = 3,
      inherit.aes = FALSE
    ) +
    geom_errorbarh(
      data = tip.rates.sum,
      aes(
        xmin = Mean - SD,
        xmax = Mean + SD,
        y = c(1.2, 2.2)
      ),
      height = 0.05,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(breaks = pretty(c(0, max(tip.rates.sum$Max)), n = 8)) +
    labs(y = "", x = parameter) +
    theme_classic() +
    theme(legend.position = "none")
}

#' Plot HiSSE model-averaged marginal ancestral state reconstruction for the trait
#'
#' @description A function to plot a (model-averaged) marginal ancestral reconstruction for the trait data.
#'
#' @param processed_hisse_recon An object produced by \code{h_process_recon}
#' @param x_label The name of the trait to be used for guide title
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
#'    processed_hisse_recon = processed_hisse,
#'    x_lab = "", discrete=FALSE, cutoff=.5)
#'
#'# change colors, your can pass the trait name to `name=` to title the colorbar
#'map_continuous + scale_color_gradient(name="", low = "#132B43", high = "#56B1F7")
#'
#'map_discrete <-
#'  h_trait_recon(
#'    processed_hisse_recon = processed_hisse,
#'    x_lab = "", discrete=TRUE, cutoff=.5)
#'
#'# change colors
#'map_discrete + scale_color_manual(name="", values = c("red", "blue"))
#'
#'@export

h_trait_recon <-
  function(processed_hisse_recon,
           x_label,
           discrete = FALSE,
           cutoff = .5,
           tree_layout = "rectangular",
           tree_direction = "right",
           time_axis_ticks = 10,
           open_angle = 10) {
    if (!tree_layout %in% c('rectangular', 'circular', 'slanted', 'fan', 'radial')) {
      stop("The selected tree layout is not supported.")
    }

    tree <- processed_hisse_recon$tree_data@phylo
    state <- processed_hisse_recon$tree_data@data$state
    agemax <- tree %>% branching.times() %>% max()

    ggg <-
      ggtree(
        tr = tree,
        layout = tree_layout,
        size = .45,
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
      d_state <- ifelse(state > cutoff, 1, 0) %>% as.factor()
      ggg <- ggg +
        aes(color = d_state) +
        scale_color_viridis_d(name = x_label, end = 0.6) +
        guides(color = guide_legend(override.aes = list(size = 4)))
    } else {
      d_state <- state
      ggg <- ggg +
        aes(color = d_state) +
        scale_color_viridis_c(name = x_label, end = 0.6)
    }

    ggg <-
      tree_flip(
        ggtree_object = ggg,
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
#' @param processed_hisse_recon An object produced by \code{h_process_recon}
#' @param parameter The diversification parameter to be mapped onto the tree. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param discrete Logical. Whether to discretize the distribution of reconstructed rates into bins
#' @param breaks A numeric vector of cut points for binning the rates. Passed internally to \code{cut}
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
#'  h_rate_recon(
#'    processed_hisse_recon = processed_hisse,
#'    parameter = "extinction", discrete=FALSE)
#'
#'# change colors, your can pass the trait name to `name=` to title the colorbar
#'map_continuous + scale_color_gradient(name="", low = "#132B43", high = "#56B1F7")
#'
#'map_discrete <-
#'  h_rate_recon(
#'    processed_hisse_recon = processed_hisse,
#'    parameter = "extinction", discrete=TRUE, breaks=c(0.3, 0.6, 1))
#'
#'# change colors
#'map_discrete + scale_color_manual(name="", values = c("red", "blue", "orange", "green"))
#'
#'@export

h_rate_recon <-
  function(processed_hisse_recon,
           parameter = "turnover",
           discrete = FALSE,
           breaks = seq(0, 1, 0.2),
           tree_layout = "rectangular",
           tree_direction = "right",
           time_axis_ticks = 10,
           open_angle = 10) {
    if (!tree_layout %in% c('rectangular', 'circular', 'slanted', 'fan', 'radial')) {
      stop("The selected tree layout is not supported.")
    }

    tree <- processed_hisse_recon$tree_data@phylo
    datas <- processed_hisse_recon$tree_data@data
    agemax <- tree %>% branching.times() %>% max()
    wanted <- as.name(parameter)

    ggg <-
      ggtree(
        tr = tree,
        layout = tree_layout,
        size = 0.45,
        open.angle = 10
      ) +
      theme(
        # legend.direction = "horizontal",
        # legend.position = "bottom",
        legend.key.size = unit(x = 0.5, units = "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank()
      )

    if (discrete) {
      param <- datas %>%
        select(!!wanted) %>%
        unlist %>%
        unname %>%
        cut(breaks = breaks)
      ggg <-
        ggg + aes(color = param) + scale_color_viridis_d(
          option = "B",
          begin = 0.25,
          end = 0.8,
          name = parameter
        ) + guides(color = guide_legend(override.aes = list(size = 4)))
    } else {
      param <- datas %>% select(!!wanted) %>% unlist %>% unname
      ggg <-
        ggg + aes(color = param) + scale_color_viridis_c(
          option = "B",
          begin = 0.25,
          end = 0.8,
          name = parameter
        )
    }

    ggg <-
      tree_flip(
        ggtree_object = ggg,
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
#'States <-
#'    expand.grid(observed=c("marine-plankton",
#'                           "marine-benthos",
#'                           "freshwater-plankton",
#'                           "freshwater-benthos"),
#'               hidden=LETTERS[1:2]) %>%
#'    mutate(Name=paste(observed, hidden, sep="_")) %>%
#'    select(Name) %>% unlist %>% unname
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
    Dim <- dim(model_fit$trans.matrix)[1]

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
      wanted <-
        colnames(model_fit$trans.matrix) %>%
        str_replace(regex("\\("), "") %>%
        str_replace(regex("\\)"), "") %>%
        expand.grid(., .) %>%
        mutate("Names" = paste("q", Var1, "_", Var2, sep = "")) %>%
        select(Names) %>% unlist %>% unname
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
      unlist %>% zapsmall(x = ., digits = 8) %>% matrix(ncol = Dim)

    Col_names <-
      colnames(model_fit$trans.matrix) %>%
      str_replace(regex("\\("), "") %>%
      str_replace(regex("\\)"), "")

    translation <- data.frame(Col_names, states_names = states)
    colnames(Rate_matrix) <-
      rownames(Rate_matrix) <- translation$states_names
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
  hidden_traits <- ifelse(hidden_traits == "A", 1, 2)
  num_hidden_traits <- length(hidden_traits)
  num_hidden_traits <- num_hidden_traits * 4
  tur <- model_fit$solution %>%
    t %>%
    data.frame %>%
    select(starts_with("turnover")) %>%
    .[, 1:num_hidden_traits] %>%
    unname %>%
    zapsmall
  eps <- model_fit$solution %>%
    t %>%
    data.frame %>%
    select(starts_with("eps")) %>%
    .[, 1:num_hidden_traits] %>%
    unname %>%
    zapsmall
  States <- expand.grid(states, hidden_traits) %>%
    mutate("Name" = paste(Var1, Var2, sep = "")) %>%
    select(Name)

  # if (transform) {
  #   spec_ext <- hisse:::ParameterTransform(x = tur, y = eps)
  #   res <- data.frame(State=States, Speciation=t(spec_ext[1:length(tur)]), Extinction=t(spec_ext[(length(tur)+1):length(spec_ext)]))
  # } else {
  #   res <- data.frame(State=States, Turnover=t(tur), Extinction_fraction=t(eps))
  # }

  res <-
    data.frame(
      State = States,
      Turnover = t(tur),
      Extinction_fraction = t(eps)
    )
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
    mat_size <- model_fit$trans.matrix %>% ncol
    states <-
      expand.grid(character_states, 1:(mat_size / 4)) %>%
      mutate(Name = paste(Var1, Var2, sep = "")) %>%
      select(Name) %>% unlist %>% unname

    trans_rates <-
      m_transition_matrix(model_fit = model_fit,
                          hidden_states = hidden_traits,
                          states = states)
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
    mutate("prob_0x" = state.00 + state.01,
           "prob_x0" = state.00 + state.10)
  colnames(tip.rates)[1] <- "id"

  nod.rates <-
    GetModelAveRates(x = muhisse_recon, type = "nodes") %>%
    as_tibble() %>%
    mutate("prob_0x" = state.00 + state.01,
           "prob_x0" = state.00 + state.10)

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
#' @param processed_muhisse_recon An object produced with \code{m_process_recon}
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
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'m_scatterplot_cp(
#'  processed_muhisse_recon = processed_muhisse,
#'  parameter = "turnover",
#'  focal_character = "prob_0x",
#'  focal_character_label = "p(mar)",
#'  second_character_label = "p(pla)",
#'  colors = viridis(n = 9)[c(5, 1, 9)],
#'  plot_as_waiting_time = TRUE) +
#'  labs(y="Net turnover\n(waiting time in millions of years)")
#'
#'@export

m_scatterplot_cp <-
  function(processed_muhisse_recon,
           parameter = "turnover",
           focal_character = c("prob_0x", "prob_x0"),
           focal_character_label,
           second_character_label,
           colors,
           plot_as_waiting_time = FALSE) {
    if (plot_as_waiting_time) {
      tip_rates <- processed_muhisse_recon$tip_rates %>%
        mutate(wanted = 1 / !!as.name(parameter))
    } else
      (tip_rates <- processed_muhisse_recon$tip_rates %>%
         mutate(wanted = !!as.name(parameter)))

    max_rate <- tip_rates %>% select(wanted) %>% unlist %>% max

    sum_tip_rates <- tip_rates %>%
      mutate(both_prob = interaction(prob_0x, prob_x0)) %>%
      group_by(both_prob) %>%
      select(both_prob, wanted) %>%
      summarise_if(is.numeric, .funs = list("MN" = mean, "SD" = sd)) %>%
      mutate(LB = MN - SD, UB = MN + SD) %>%
      mutate(prob_0x = str_extract(both_prob, regex("^\\d+")) %>% as.numeric()) %>%
      mutate(prob_x0 = str_extract(both_prob, regex("\\d+$")) %>% as.numeric())

    if (focal_character == "prob_0x") {
      sum_tip_rates <-
        mutate(sum_tip_rates, focal_character = factor(prob_0x))
      print(sum_tip_rates)
      sss <-
        ggplot(tip_rates,
               aes(
                 x = factor(!!as.name(focal_character)),
                 y = wanted,
                 color = prob_0x,
                 color2 = prob_x0
               ))
      nudgex <- c(-0.3, -0.3, 0.3, 0.3)
    }

    if (focal_character == "prob_x0") {
      sum_tip_rates <-
        mutate(sum_tip_rates, focal_character = factor(prob_x0))
      print(sum_tip_rates)
      sss <-
        ggplot(tip_rates,
               aes(
                 x = factor(!!as.name(focal_character)),
                 y = wanted,
                 color = prob_x0,
                 color2 = prob_0x
               ))
      nudgex <- c(-0.3, 0.3, -0.3, 0.3)
    }

    sss <- sss +
      geom_point(alpha = .7,
                 size = 1.25,
                 position = position_jitter(width = .15)) +
      geom_errorbar(
        data = sum_tip_rates,
        aes(
          x = focal_character,
          y = MN,
          ymin = LB,
          ymax = UB
        ),
        position = position_nudge(x = nudgex, y = 0),
        inherit.aes = TRUE,
        width = .04
      ) +
      geom_point(
        data = sum_tip_rates,
        aes(x = focal_character, y = MN),
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
      labs(x = focal_character_label, y = parameter)
    return(sss)
  }


#' Plot diversification rates estimated by a MuHiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a ridgeline of (model-averaged) diversification rates in the alternative states.
#'
#' @param processed_muhisse_recon An object produced by \code{m_process_recon}
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
#'  processed_muhisse_recon = processed_muhisse,
#'  parameter = "extinction",
#'  line_colors = viridis(n = 4, option=1))
#'
#'@export

m_ridgelines <- function(processed_muhisse_recon,
                         states_names = c("00", "01", "10", "11"),
                         parameter = "turnover",
                         plot_as_waiting_time = FALSE,
                         fill_colors = rep(NA, 4),
                         line_colors = viridis(n = 4)) {
  message(
    "Recoding and renaming character states. The elements 1:4 of the vector `character_states_names` are assumed to match the states 00, 01, 10, 11.\n\n"
  )

  ss <-
    m_prep_df(
      processed_muhisse_recon = processed_muhisse_recon,
      states_names = states_names,
      parameter = parameter
    )

  message("Summarising grouped by character state\n\n")
  wanted <- as.name(parameter)

  summ <- . %>% summarise_at(.vars = vars(wanted),
                             .funs = list(
                               "Mean" = mean,
                               "Median" = median,
                               "SD" = sd
                             ))

  if (plot_as_waiting_time) {
    ss <- mutate(ss, wanted = 1 / !!wanted)
    ss.sum <- ss %>% group_by(four_state) %>%
      summ(.)
  } else {
    ss <- mutate(ss, wanted = !!wanted)
    ss.sum <- ss %>% group_by(four_state) %>%
      summ(.)
  }
  max_rate <-
    ss %>% select(wanted) %>% top_n(1, wt = wanted) %>% unlist %>% unname
  print(ss.sum)

  message("\nPlotting\n\n")

  ggplot(data = ss,
         aes(
           x = wanted,
           y = four_state,
           fill = four_state,
           colour = four_state
         )) +
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
      position = position_nudge(y = rep(-0.35, 4)),
      aes(
        xmin = Mean - SD,
        xmax = Mean + SD,
        y = four_state,
        colour = four_state
      ),
      height = 0.05,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = ss.sum,
      pch = 21,
      position = position_nudge(y = rep(-0.3, 4)),
      aes(y = four_state,
          x = Mean,
          colour = four_state),
      size = 3,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(breaks = pretty(x = c(0, max_rate), n = 10)) +
    scale_fill_manual(values = fill_colors) +
    scale_colour_manual(values = line_colors) +
    labs(y = "", x = parameter) +
    theme_classic() +
    theme(legend.position = "none")
}

#' Plot diversification rates estimated by a MuHiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a jittered scatterplot of (model-averaged) diversification rates and summary statistics in the alternative states.
#'
#' @param processed_muhisse_recon An object produced by \code{m_process_recon}
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
#'  processed_muhisse_recon = processed_muhisse,
#'  parameter = "extinction",
#'  colors = viridis(n = 4, option=1, end=.8))
#'
#'@export

m_scatterplot <-
  function(processed_muhisse_recon,
           states_names = c("00", "01", "10", "11"),
           parameter = "turnover",
           colors = viridis(n = 4, option = 2, end = 0.7),
           plot_as_waiting_time = FALSE) {
    message(
      "Recoding and renaming character states. The elements 1:4 of the vector `character_states_names` are assumed to match the states 00, 01, 10, 11.\n"
    )

    ss <-
      m_prep_df(
        processed_muhisse_recon = processed_muhisse_recon,
        states_names = states_names,
        parameter = parameter
      )

    message("Summarising grouped by character state\n")
    wanted <- as.name(parameter)

    summ <- . %>% summarise_at(.vars = vars(wanted),
                               .funs = list(
                                 Mean = mean,
                                 Median = median,
                                 SD = sd
                               )) %>%
      mutate(wanted = Mean) # to bypass warning in geom_errorbar

    if (plot_as_waiting_time) {
      ss <- mutate(ss, wanted = 1 / !!wanted)
      ss.sum <- ss %>% group_by(four_state) %>%
        summ(.)
    } else {
      ss <- mutate(ss, wanted = !!wanted)
      ss.sum <- ss %>% group_by(four_state) %>%
        summ(.)
    }
    max_rate <-
      ss %>% select(wanted) %>% top_n(1, wt = wanted) %>% unlist %>% unname
    print(select(ss.sum,-wanted))

    message("\nPlotting\n")

    pl <-
      ggplot(data = ss,
             aes(x = four_state,
                 y = wanted,
                 fill = four_state)) +
      geom_point(
        pch = 21,
        position = position_jitter(0.2),
        size = 2,
        stroke = .3
      ) +
      geom_errorbar(
        inherit.aes = FALSE,
        data = ss.sum,
        aes(
          x = four_state,
          colour = four_state,
          ymin = Mean - SD,
          ymax = Mean + SD
        ),
        width = .05,
        position = position_nudge(x = .35)
      ) +
      geom_point(
        inherit.aes = FALSE,
        data = ss.sum,
        aes(x = four_state, y = Mean, colour = four_state),
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
      theme(legend.position = "top")
    return(pl)
  }

#' Plot diversification rates estimated by a MuHiSSE model with means and standard deviations across tips
#'
#' @description A function to plot a stacked dotplot of (model-averaged) diversification rates and summary statistics in the alternative states.
#'
#' @param processed_muhisse_recon An object produced by \code{m_process_recon}
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
#'  processed_muhisse_recon = processed_muhisse,
#'  parameter = "extinction",
#'  colors = viridis(n = 4, option=1, end=.8))
#'
#'@export

m_dotplot <-
  function(processed_muhisse_recon,
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
        processed_muhisse_recon = processed_muhisse_recon,
        states_names = states_names,
        parameter = parameter
      )

    message("Summarising grouped by character state\n\n")
    wanted <- as.name(parameter)

    summ <- . %>% summarise_at(.vars = vars(wanted),
                               .funs = list(
                                 Mean = mean,
                                 Median = median,
                                 SD = sd
                               ))

    if (plot_as_waiting_time) {
      ss <- mutate(ss, wanted = 1 / !!wanted)
      ss.sum <- ss %>% group_by(four_state) %>%
        summ(.)
    } else {
      ss <- mutate(ss, wanted = !!wanted)
      ss.sum <- ss %>% group_by(four_state) %>%
        summ(.)
    }
    max_rate <-
      ss %>% select(wanted) %>% top_n(1, wt = wanted) %>% unlist %>% unname
    # print(ss)
    print(ss.sum)

    cat ("\nPlotting\n\n")

    pl <- ggplot(data = ss,
                 aes(x = four_state,
                     y = wanted,
                     fill = four_state)) +
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
          x = four_state,
          y = Mean,
          colour = four_state,
          ymin = Mean - SD,
          ymax = Mean + SD
        ),
        width = .05,
        position = position_nudge(x = -.1)
      ) +
      geom_point(
        inherit.aes = FALSE,
        data = ss.sum,
        aes(x = four_state, y = Mean, colour = four_state),
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
      theme(legend.position = "none")
    return(pl)
  }

#' Plot MuHiSSE model-averaged marginal ancestral state reconstruction for the trait with a two-dimensional colorplane
#'
#' @description A function to plot a MuHiSSE (model-averaged) marginal ancestral reconstruction for the trait data.
#'
#' @param processed_muhisse_recon An object produced by \code{m_process_recon}
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
#'data("diatoms")
#'processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
#'
#'m_trait_recon_cp(
#'  processed_muhisse_recon = processed_muhisse,
#'  tree_layout = "fan",
#'  focal_character = "prob_0x",
#'  focal_character_label = "p(marine)",
#'  second_character_label = "p(plankton)",
#'  colors = viridis(9)[c(5,1,9)]
#'  )
#'
#'@export

m_trait_recon_cp <-
  function(processed_muhisse_recon,
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

    tree <- processed_muhisse_recon$tree_data@phylo
    datas <- processed_muhisse_recon$tree_data@data
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
#' @param processed_muhisse_recon An object produced by \code{m_process_recon}
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
#'
#'# eight categories after binning with a cutoff of 0.2
#'m_trait_recon(
#'  processed_muhisse_recon = processed_muhisse,
#'  cutoff = c(.2, .2),
#'  states_of_first_character = c("marine", "freshwater"),
#'  states_of_second_character = c("plankton", "benthos"),
#'  tree_layout = "radial",
#'  colors = viridis(7,option = "E", direction = 1, end=.9))
#'
#'# three of these eight have < 3 nodes, so we could try to avoid plotting some of them
#'# adjust the cutoff for the second variable
#'
#'m_trait_recon(
#'  processed_muhisse_recon = processed_muhisse,
#'  cutoff = c(.2, .3),
#'  states_of_first_character = c("marine", "freshwater"),
#'  states_of_second_character = c("plankton", "benthos"),
#'  tree_layout = "radial",
#'  colors = viridis(7,option = "E", direction = 1, end=.9))
#'
#'# ignoring uncertainty
#'m_trait_recon(
#'  processed_muhisse_recon = processed_muhisse,
#'  cutoff = c(.5, .5),
#'  states_of_first_character = c("marine", "freshwater"),
#'  states_of_second_character = c("plankton", "benthos"),
#'  tree_layout = "radial",
#'  colors = viridis(7,option = "E", direction = 1, end=.9)[c(1,3,5,7)])
#'
#'@export

m_trait_recon <-
  function(processed_muhisse_recon,
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

    tree <- processed_muhisse_recon$tree_data@phylo
    agemax <- tree %>% branching.times() %>% max()

    ss <- processed_muhisse_recon$tree_data@data %>%
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
      mutate(four_state = paste(prob_0x_named, prob_x0_named, sep = "-"))

    message("Categories after discretizing with the provided cutoff:\n")
    ss.cnt <- ss %>%
      ungroup %>%
      group_by(four_state) %>%
      add_tally() %>%
      select(four_state, n) %>%
      distinct

    print(ss.cnt)

    nstat <- ss %>% select(four_state) %>% distinct %>% nrow

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
#' @param processed_muhisse_recon An object produced by \code{h_process_recon}
#' @param parameter The diversification parameter to be mapped onto the tree. Possible options are turnover, extinct.frac, net.div, speciation, extinction
#' @param discrete Logical. Whether to discretize the distribution of reconstructed rates into bins
#' @param breaks A numeric vector of cut points for binning the rates. Passed internally to \code{cut}. The function checks whether max(rate) is in this vector and adds it if not.
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
#'map_continuous <-
#'  m_rate_recon(
#'    processed_muhisse_recon = processed_muhisse,
#'    parameter = "extinction", discrete=FALSE)
#'
#'# change colors, your can pass the trait name to `name=` to title the colorbar
#'map_continuous + scale_color_gradient(name="", low = "#132B43", high = "#56B1F7")
#'
#'map_discrete <-
#'  m_rate_recon(
#'    processed_muhisse_recon = processed_muhisse,
#'    parameter = "extinction", discrete=TRUE, breaks=c(0.3, 0.6, 1))
#'
#'# change colors
#'map_discrete + scale_color_manual(name="", values = c("red", "blue", "orange", "green"))
#'
#'@export

m_rate_recon <-
  function(processed_muhisse_recon,
           parameter = "turnover",
           discrete = FALSE,
           breaks = seq(0, 1, 0.2),
           tree_layout = "rectangular",
           tree_direction = "right",
           time_axis_ticks = 10,
           open_angle = 10) {
    if (!tree_layout %in% c('rectangular', 'circular', 'slanted', 'fan', 'radial')) {
      stop("The selected tree layout is not supported.")
    }

    tree <- processed_muhisse_recon$tree_data@phylo
    datas <- processed_muhisse_recon$tree_data@data
    agemax <- tree %>% branching.times() %>% max()
    wanted <- as.name(parameter)

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
      max_rate <- datas %>% select(!!wanted) %>% unlist %>% unname %>% max
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
        select(!!wanted) %>%
        unlist %>%
        unname %>%
        cut(breaks = breaks)
      ggg <-
        ggg + aes(color = param) + scale_color_viridis_d(
          option = "B",
          begin = 0.25,
          end = 0.8,
          name = parameter
        ) + guides(color = guide_legend(override.aes = list(size = 4)))
    } else {
      param <- datas %>% select(!!wanted) %>% unlist %>% unname
      ggg <-
        ggg + aes(color = param) + scale_color_viridis_c(
          option = "B",
          begin = 0.25,
          end = 0.8,
          name = parameter
        )
    }
    ggg <-
      tree_flip(
        ggtree_object = ggg,
        tree_layout = tree_layout,
        tree_direction = tree_direction,
        time_axis_ticks = time_axis_ticks,
        agemax = agemax
      )

    return(ggg + theme(plot.margin = unit(rep(.1, 4), "in")))
  }


##### --- Utility functions ------------------------- #####

tree_flip <- function(ggtree_object,
                      tree_layout,
                      tree_direction,
                      time_axis_ticks,
                      agemax) {
  if (tree_layout %in% c("rectangular", "slanted")) {
    if (tree_direction == "up") {
      ggtree_object <- ggtree_object + coord_flip() +
        theme(
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.text.y = element_text()
        ) +
        scale_x_continuous(
          expand = c(0, 0.01),
          breaks = pretty(c(0, agemax), n = time_axis_ticks),
          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))
        ) +
        labs(x = "Time (Ma)")
    }
    if (tree_direction == "down") {
      ggtree_object <- ggtree_object + coord_flip() +
        theme(
          axis.line.y = element_line(),
          axis.ticks.y = element_line(),
          axis.text.y = element_text()
        ) +
        scale_x_continuous(
          trans = "reverse",
          expand = c(0, 0.01),
          breaks = pretty(c(0, agemax), n = time_axis_ticks),
          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))
        ) +
        labs(x = "Time (Ma)")
    }

    if (tree_direction == "left") {
      ggtree_object <- ggtree_object +
        theme(
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text()
        ) +
        scale_x_continuous(
          trans = "reverse",
          expand = c(0, 0.01),
          breaks = pretty(c(0, agemax), n = time_axis_ticks),
          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))
        ) +
        labs(x = "Time (Ma)")
    }

    if (tree_direction == "right") {
      ggtree_object <- ggtree_object +
        theme(
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text()
        ) +
        scale_x_continuous(
          expand = c(0, 0.01),
          breaks = pretty(c(0, agemax), n = time_axis_ticks),
          labels = rev(pretty(c(0, agemax), n = time_axis_ticks))
        ) +
        labs(x = "Time (Ma)")
    }
  }

  if (tree_layout %in% c("circular", "fan", "radial")) {
    maxx <- ggtree_object$data %>%
      top_n(n = 1, wt = x) %>%
      select(x) %>%
      unlist %>%
      unname %>%
      unique
    maxx <- round(maxx, 1)
    ntip <- ggtree_object$data %>% filter(isTip==TRUE) %>% nrow()
    ntip <- ntip+10
    pretty_points <-
      maxx - c(maxx, pretty(c(maxx:0), n = time_axis_ticks))

    pp <- tibble(x = rev(pretty_points), y = 0) %>%
      filter(x <= maxx, x > 0) %>%
      mutate(label = rev(x) - min(x))

    ggtree_object <- ggtree_object +
      geom_vline(
        data = pp,
        aes(xintercept = x),
        size = .2,
        color = "darkgrey"
      ) +
      geom_text(
        data = pp,
        aes(
          x = x + 0.1,
          y = ntip + 2,
          label = label
        ),
        size = 2,
        inherit.aes = FALSE
      )
  }
  return(ggtree_object)
}

m_prep_df <-
  function(processed_muhisse_recon,
           states_names,
           parameter) {
    ss <- processed_muhisse_recon$tip_rates %>%
      mutate(what_state = paste(state.00, state.01, state.10, state.11, sep = "")) %>%
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
