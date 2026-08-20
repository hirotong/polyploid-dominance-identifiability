## Internal null-coalescing operator. Not exported: base R only gained `%||%`
## in 4.4, but DESCRIPTION declares Depends: R (>= 4.0), so polyDx defines its
## own. Used by later polyDx code; kept here so it is defined before anything
## else in the package loads.
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Read a polyDx run-time option from the environment
#'
#' The workflow is driven on the cluster by environment variables
#' (\code{PWF_NREF}, \code{PWF_FRAC_NITER}, ...). \code{pwf_opt} reads one and
#' coerces it to the type of \code{default}, falling back to \code{default} when
#' the variable is unset, empty, or fails to coerce.
#'
#' When \code{default} is logical, the coercion rule is \code{raw != "0"}:
#' the literal string \code{"0"} coerces to \code{FALSE}, and every other
#' non-empty value -- including \code{"FALSE"} -- coerces to \code{TRUE}.
#' This matches the historical driver's convention (e.g.
#' \code{PWF_MARKER_QC=FALSE} still \emph{enables} marker QC; only
#' \code{PWF_MARKER_QC=0} disables it) and is preserved for parity with
#' \code{main}, not a bug.
#'
#' @param name environment variable name.
#' @param default fallback value; its type drives the coercion.
#' @return the coerced option value.
#' @examples
#' pwf_opt("PWF_NREF", 30L)
#' @export
pwf_opt <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  out <- if (is.integer(default)) suppressWarnings(as.integer(raw))
         else if (is.numeric(default)) suppressWarnings(as.numeric(raw))
         else if (is.logical(default)) raw != "0"
         else raw
  if (length(out) != 1L || (!is.character(out) && !is.finite(out))) default else out
}

#' Default heritability for panel-conditioned simulations
#'
#' All panel-conditioned controls simulate at this broad-sense heritability
#' unless a panel-specific value is supplied.
#'
#' @format numeric scalar, 0.6.
#' @examples
#' PWF_H2_DEFAULT
#' @export
PWF_H2_DEFAULT <- 0.6

#' Decision thresholds for the diagnostic scorecard
#'
#' The gate cut-offs used by \code{\link{workflow_scorecard}}. These are the
#' values that produced the manuscript's scorecard; override them by passing a
#' modified list.
#'
#' @return a named list with \code{struct_large}, \code{recover_rate},
#'   \code{ci_wide}, \code{null_margin} and \code{pred_gain}.
#' @examples
#' workflow_thresholds()
#' @export
workflow_thresholds <- function() {
  list(
    struct_large = 0.10,  # |S_struct_real| at or above this => structure-sensitive
    recover_rate = 0.80,  # D2 spike recovery rate at or above this => recoverable
    ci_wide      = 0.25,  # D2 fraction CrI width at or above this => too wide
    null_margin  = 0.05,  # min delta_null to count as exceeding the panel's null
    pred_gain    = 0.01   # material mean paired CV gain of A+D over A
  )
}
