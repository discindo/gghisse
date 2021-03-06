% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilhisse.R
\name{h_ridgelines}
\alias{h_ridgelines}
\title{Plot diversification rates estimated by a HiSSE model with means and standard deviations across tips}
\usage{
h_ridgelines(
  processed_recon,
  parameter = "turnover",
  states_names = c("Marine", "Freshwater"),
  line_colors = c("black", "black"),
  fill_colors = c("yellow", "red"),
  plot_as_waiting_time = FALSE
)
}
\arguments{
\item{processed_recon}{An object produced by \code{h_process_recon}}

\item{parameter}{The diversification parameter to be plotted on the y axis. Possible options are turnover, extinct.frac, net.div, speciation, extinction}

\item{states_names}{The names for character states}

\item{line_colors}{Colors for the lines and points}

\item{fill_colors}{Colors for the density polygons}

\item{plot_as_waiting_time}{Whether to plot the rates (FALSE, default) or their inverse (waiting times)}
}
\value{
A ridgeline plot of tip-associated rates (possibly model averaged).
}
\description{
A function to plot a ridgeline of (model-averaged) diversification rates in the alternative states. We can change the rate plotted on the y axis and modify the label for the x-axis (your binary 0/1 trait). Further modifications are straightforward with ggplot (see examples).
}
\examples{

data("diatoms")
processed_hisse <- h_process_recon(hisse_recon=diatoms$cid4_recon)

h_ridgelines(
 processed_recon = processed_hisse,
 states_names = c("Plankton", "Benthos"),
 parameter = "extinction",
 line_colors = c("black", "black"),
 fill_colors = c("orange", "violet"))

}
