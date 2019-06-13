#' Some hisse objects, either from a diatom dataset or simulated.
#'
#'
#' @format A list with 5 elements
#' \describe{
#'   \item{cid4_recon}{hisse.states object. A marginal ancestral state reconstruction of diatom habitat (plankton-benthos) under hisse's CID4 model}
#'   \item{cid8}{muhisse.fit object. A CID8 model fit object for data on diatom salinity (marine-freshwater) and habitat (plankton-benthos)}
#'   \item{muhisse}{muhisse.fit object. A MuHiSSE model fit object for data on diatom salinity (marine-freshwater) and habitat (plankton-benthos)}
#'   \item{muhisse_recon}{muhisse.states object. A marginal ancestral state reconstruction of diatom salinity (marine-freshwater) and habitat (plankton-benthos) under hisse's MuHiSSE model}
#'   \item{musse}{muhisse.fit object. A MuSSE model fit object for data on diatom salinity (marine-freshwater) and habitat (plankton-benthos)}
#'   \item{3state_musse_fit}{muhisse.fit object. A MuSSE model for simulated data with three states}
#'   \item{3state_musse_recon}{muhisse.states object. A MuSSE marginal ancestral state reconstruction for simulated data with three states}
#' }
#' @references \url{https://www.biorxiv.org/content/10.1101/406165v2.abstract}
"diatoms"
