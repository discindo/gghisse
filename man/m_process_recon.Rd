% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilhisse.R
\name{m_process_recon}
\alias{m_process_recon}
\title{Prepare MuHiSSE marginal reconstruction for plotting}
\usage{
m_process_recon(muhisse_recon)
}
\arguments{
\item{muhisse_recon}{An object produced by \code{hisse::MarginReconMuHiSSE}, or a list of such objects over which we can model-average using \code{hisse::GetModelAveRates}}
}
\description{
A function to process the output from \code{hisse::MarginReconMuHiSSE}. First we get model averaged rates for tips and nodes of the phylogeny and probabilities from the ancestral reconstruction of character states. Then convert the tree from \code{phylo} format into a \code{treedata} class and assign the data.frame of mode-averaged rates to the \code{@data} slot. We output the node and tip rates tables and the reformated tree object for downstream plotting with \code{ggtree} and \code{ggplot}.
}
\details{
Internally, the probabilities for individual states (00, 01, 10, 11) are converted into \code{prob_0x} giving the probability for the first character being 0, and \code{prob_x0} giving the probability of the second character being 0. This way we can plot in terms of a 'focal character', i.e., we can plot the rate on the y axis and the probability of being 0 for the 'focal character' on the x axis while color-coding by the states of both characters.

#'@return A list with three components: \itemize{
 \item \code{tree_data} the phylogeny with associated data object with reconstructed ancestral, model-averaged node and tip states and rates
 \item \code{tip_rates/node_rates} a tibble of reconstructed ancestral, model-averaged tip/node states and rates.
 }
}
\examples{

data("diatoms")
processed_muhisse <- m_process_recon(muhisse_recon=diatoms$muhisse_recon)
processed_muhisse$tip_rates
processed_muhisse$node_rates
processed_muhisse$tree_data

}
