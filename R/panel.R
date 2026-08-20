#' Parse a ploidy specification
#'
#' Accepts a fixed even ploidy as a string (\code{"4"}, \code{"6"}) or a
#' per-marker variable-ploidy spec \code{"var:<marker_ploidy.rds>:<maxdeg>"}
#' (e.g. sugarcane), where the RDS holds one integer ploidy per marker.
#'
#' @param arg the specification string.
#' @param nmarkers number of markers, used to validate a variable-ploidy vector.
#' @return a list with \code{mode} (\code{"fixed"} or \code{"variable"}),
#'   \code{ploidy} or \code{mploidy}, \code{maxdeg}, \code{ndeg} and \code{label}.
#' @examples
#' parse_ploidy("4", 100L)
#' @export
parse_ploidy <- function(arg, nmarkers) {
  if (grepl("^var:", arg)) {
    parts <- strsplit(arg, ":", fixed = TRUE)[[1]]
    if (length(parts) < 2 || !nzchar(parts[2])) {
      stop("variable ploidy syntax is var:<marker_ploidy.rds>:<maxdeg>", call. = FALSE)
    }
    mp <- as.integer(readRDS(parts[2]))
    if (length(mp) != nmarkers) stop("marker ploidy vector length must equal ncol(dosage)", call. = FALSE)
    maxdeg <- if (length(parts) >= 3 && nzchar(parts[3])) as.integer(parts[3]) else min(mp)
    maxdeg <- min(maxdeg, min(mp))
    return(list(mode = "variable", mploidy = mp, maxdeg = maxdeg, ndeg = maxdeg,
                label = paste0("var<=", maxdeg)))
  }
  pl <- as.integer(arg)
  if (!is.finite(pl) || pl < 2) stop("fixed ploidy must be an integer >= 2", call. = FALSE)
  list(mode = "fixed", ploidy = pl, maxdeg = pl, ndeg = pl, label = as.character(pl))
}

#' Subset a ploidy specification to a retained marker set
#'
#' @param pcfg a \code{\link{parse_ploidy}} result.
#' @param keep logical or integer index over the original markers.
#' @return \code{pcfg} with \code{mploidy} subset (variable mode) or unchanged.
#' @examples
#' subset_ploidy(parse_ploidy("4", 3L), c(TRUE, FALSE, TRUE))
#' @export
subset_ploidy <- function(pcfg, keep) {
  if (pcfg$mode == "variable") pcfg$mploidy <- pcfg$mploidy[keep]
  pcfg
}

#' Marker QC: call rate, MAF, polymorphism, then mean imputation
#'
#' Drops markers below the call-rate and MAF thresholds and monomorphic markers,
#' mean-imputes the remaining missing cells, rounds to integer dosage classes and
#' clamps to \code{0..ploidy}. Returns a one-row QC audit alongside the matrix.
#'
#' @param M individuals x markers dosage matrix (\code{NA} allowed).
#' @param pcfg a \code{\link{parse_ploidy}} result.
#' @param maf_min minimum minor-allele frequency, a single finite value in
#'   \code{[0, 0.5]}; defaults to the \code{PWF_MAF_MIN} environment variable
#'   via \code{\link{pwf_opt}}, or 0.05. Validated at call time.
#' @param call_min minimum call rate (strict \code{>}), a single finite value
#'   in \code{[0, 1]}; defaults to the \code{PWF_CALLRATE_MIN} environment
#'   variable via \code{\link{pwf_opt}}, or 0.8. Validated at call time.
#' @param enabled if \code{FALSE}, keep every marker (imputation still applies);
#'   defaults to the \code{PWF_MARKER_QC} environment variable via
#'   \code{\link{pwf_opt}} (\code{"0"} disables, anything else including unset
#'   enables).
#' @return a list with \code{M}, \code{pcfg} and a one-row \code{qc} data frame.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(400, 4, 0.4), nrow = 20)
#' filter_markers(M, parse_ploidy("4", ncol(M)))$qc$markers_post_qc
#' @export
filter_markers <- function(M, pcfg,
                           maf_min = pwf_opt("PWF_MAF_MIN", 0.05),
                           call_min = pwf_opt("PWF_CALLRATE_MIN", 0.8),
                           enabled = pwf_opt("PWF_MARKER_QC", TRUE)) {
  if (!is.numeric(maf_min) || length(maf_min) != 1L || !is.finite(maf_min) ||
      maf_min < 0 || maf_min > 0.5) {
    stop("maf_min must be a single finite value in [0, 0.5]", call. = FALSE)
  }
  if (!is.numeric(call_min) || length(call_min) != 1L || !is.finite(call_min) ||
      call_min < 0 || call_min > 1) {
    stop("call_min must be a single finite value in [0, 1]", call. = FALSE)
  }
  n0 <- nrow(M); m0 <- ncol(M)
  missing_pre_qc <- sum(is.na(M))
  callrate <- colMeans(!is.na(M))
  if (pcfg$mode == "fixed") {
    denom <- rep(pcfg$ploidy, m0)
  } else {
    denom <- pcfg$mploidy
  }
  af <- colMeans(M, na.rm = TRUE) / denom
  maf <- pmin(af, 1 - af)
  polymorphic <- apply(M, 2, function(z) length(unique(stats::na.omit(round(z)))) > 1)
  keep <- is.finite(callrate) & callrate > call_min & is.finite(maf) &
    maf >= maf_min & polymorphic
  if (!enabled) keep <- rep(TRUE, m0)
  if (!any(keep)) stop("no markers remain after marker QC", call. = FALSE)

  M <- M[, keep, drop = FALSE]
  pcfg <- subset_ploidy(pcfg, keep)
  missing_after_filter <- sum(is.na(M))
  imputed_cells <- 0L
  if (missing_after_filter > 0) {
    for (j in which(colSums(is.na(M)) > 0)) {
      mu <- mean(M[, j], na.rm = TRUE)
      M[is.na(M[, j]), j] <- mu
    }
    imputed_cells <- missing_after_filter
  }
  M <- round(M)
  if (pcfg$mode == "fixed") {
    M[M < 0] <- 0
    M[M > pcfg$ploidy] <- pcfg$ploidy
  } else {
    MP <- matrix(pcfg$mploidy, nrow(M), ncol(M), byrow = TRUE)
    M[M < 0] <- 0
    too_high <- M > MP
    M[too_high] <- MP[too_high]
  }
  storage.mode(M) <- "double"

  qc <- data.frame(
    marker_qc_enabled = enabled,
    maf_min = maf_min,
    callrate_min = call_min,
    n_after_alignment = n0,
    markers_pre_qc = m0,
    missing_cells_pre_qc = missing_pre_qc,
    markers_callrate_pass = sum(callrate > call_min, na.rm = TRUE),
    markers_polymorphic = sum(polymorphic, na.rm = TRUE),
    markers_maf_pass = sum(maf >= maf_min, na.rm = TRUE),
    markers_post_qc = ncol(M),
    markers_dropped_qc = m0 - ncol(M),
    missing_cells_after_filter = missing_after_filter,
    imputed_cells = imputed_cells,
    min_maf_pre_qc = suppressWarnings(min(maf, na.rm = TRUE)),
    median_maf_pre_qc = suppressWarnings(stats::median(maf, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  list(M = M, pcfg = pcfg, qc = qc)
}

#' Build a panel from in-memory dosage, phenotype and group objects
#'
#' Aligns the dosage matrix and phenotype table on shared IDs, validates the
#' dosage range against the ploidy spec, applies marker QC
#' (\code{\link{filter_markers}}) and selects usable traits. This is the shared
#' body behind \code{\link{load_panel}}; call it directly when the dosage /
#' phenotype / group objects are already in memory rather than on disk.
#'
#' @param M individuals x markers dosage matrix; row names (or the default
#'   \code{id1, id2, ...}) are used as individual IDs.
#' @param phe a data frame whose first column is the ID, or \code{NULL} to skip
#'   phenotype alignment and trait selection.
#' @param ploidy_arg a \code{\link{parse_ploidy}} specification string.
#' @param group a character/factor vector positionally aligned to \code{M}'s
#'   input row order (or \code{NULL}) giving a grouping label per individual;
#'   \code{NA} entries are allowed and are re-labelled \code{"NA_grp"}.
#' @param traits \code{"auto"} or a comma-separated trait list.
#' @param max_traits maximum number of traits to keep.
#' @param min_n minimum number of overlapping individuals / finite trait values.
#' @return a list with \code{M}, \code{phe}, \code{traits}, \code{group},
#'   \code{pcfg}, \code{n}, \code{markers} and \code{qc}.
#' @seealso \code{\link{load_panel}} for reading these objects from RDS files.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(30 * 20, 4, 0.4), 30, 20)
#' rownames(M) <- paste0("id", seq_len(30))
#' phe <- data.frame(id = rownames(M), yield = rnorm(30), stringsAsFactors = FALSE)
#' P <- panel_from_matrix(M, phe, "4", min_n = 20L)
#' P$traits
#' @export
panel_from_matrix <- function(M, phe, ploidy_arg, group = NULL, traits = "auto",
                              max_traits = 8L, min_n = 40L) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  if (is.null(rownames(M))) rownames(M) <- paste0("id", seq_len(nrow(M)))
  raw_ids <- rownames(M)
  raw_n <- nrow(M); raw_markers <- ncol(M)
  pcfg <- parse_ploidy(ploidy_arg, ncol(M))

  # dosage range sanity
  if (pcfg$mode == "fixed") {
    if (any(M < 0 | M > pcfg$ploidy, na.rm = TRUE)) {
      stop(sprintf("dosage values outside 0..%d for fixed-ploidy panel", pcfg$ploidy), call. = FALSE)
    }
  } else {
    mx <- apply(M, 2, max, na.rm = TRUE)
    if (any(mx > pcfg$mploidy)) stop("some dosage values exceed their marker-specific ploidy", call. = FALSE)
  }

  traits_out <- character(0)
  group_in <- group
  if (!is.null(phe)) {
    if (is.matrix(phe)) phe <- as.data.frame(phe)
    idcol <- colnames(phe)[1]; rownames(phe) <- as.character(phe[[idcol]])
    common <- intersect(rownames(M), rownames(phe))
    if (length(common) < min_n) stop(sprintf("only %d IDs overlap dosage<->pheno (<%d)", length(common), min_n), call. = FALSE)
    if (!is.null(group_in)) group_in <- group_in[match(common, raw_ids)]
    M <- M[common, , drop = FALSE]; phe <- phe[common, , drop = FALSE]
    fq <- filter_markers(M, pcfg)
    M <- fq$M; pcfg <- fq$pcfg; qc <- fq$qc
    cand <- if (identical(traits, "auto")) setdiff(colnames(phe), idcol) else strsplit(traits, ",")[[1]]
    drop <- grepl("group|gen[ae]tic|subpop|^pop$|cluster|^id$|clone|geno", cand, ignore.case = TRUE)
    cand <- cand[!drop]
    cand <- cand[vapply(cand, function(t) {
      v <- suppressWarnings(as.numeric(phe[[t]])); is.finite(stats::sd(v, na.rm = TRUE)) && stats::sd(v, na.rm = TRUE) > 0 &&
        sum(is.finite(v)) >= min_n
    }, logical(1))]
    traits_out <- if (length(cand) > max_traits) cand[seq_len(max_traits)] else cand
  } else {
    fq <- filter_markers(M, pcfg)
    M <- fq$M; pcfg <- fq$pcfg; qc <- fq$qc
  }

  group_out <- NULL
  if (!is.null(group_in)) {
    g <- as.character(group_in)
    if (sum(!is.na(g)) > 0 && length(unique(g[!is.na(g)])) >= 2) {
      g[is.na(g)] <- "NA_grp"; group_out <- g
    }
  }

  qc$raw_n <- raw_n
  qc$raw_markers <- raw_markers
  qc$markers <- ncol(M)
  list(M = M, phe = phe, traits = traits_out, group = group_out, pcfg = pcfg,
       n = nrow(M), markers = ncol(M), qc = qc)
}

#' Load a panel from RDS files
#'
#' Reads a dosage matrix, an optional phenotype table and an optional group
#' table, aligns them on shared IDs, validates the dosage range against the
#' ploidy spec, applies marker QC and selects usable traits.
#'
#' @param dosage_rds path to an RDS holding an individuals x markers matrix with
#'   individual IDs as row names.
#' @param pheno_rds path to an RDS holding a data frame whose first column is the
#'   ID, or \code{"none"}.
#' @param ploidy_arg a \code{\link{parse_ploidy}} specification string.
#' @param group_rds path to an RDS whose first two columns are ID and group, or
#'   \code{"none"}.
#' @param traits_arg \code{"auto"} or a comma-separated trait list.
#' @param max_traits maximum number of traits to keep.
#' @param min_n minimum number of overlapping individuals / finite trait values.
#' @return a list with \code{M}, \code{phe}, \code{traits}, \code{group},
#'   \code{pcfg}, \code{n}, \code{markers} and \code{qc}.
#' @seealso \code{\link{panel_from_matrix}} for in-memory inputs.
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(30 * 20, 4, 0.4), 30, 20)
#' rownames(M) <- paste0("id", seq_len(30))
#' phe <- data.frame(id = rownames(M), yield = rnorm(30), stringsAsFactors = FALSE)
#' d <- tempfile(); dir.create(d)
#' saveRDS(M, file.path(d, "dosage.rds"))
#' saveRDS(phe, file.path(d, "pheno.rds"))
#' P <- load_panel(file.path(d, "dosage.rds"), file.path(d, "pheno.rds"), "4", min_n = 20L)
#' P$traits
#' unlink(d, recursive = TRUE)
#' @export
load_panel <- function(dosage_rds, pheno_rds, ploidy_arg, group_rds = "none",
                       traits_arg = "auto", max_traits = 8L, min_n = 40L) {
  M <- readRDS(dosage_rds)
  if (is.data.frame(M)) M <- as.matrix(M)
  storage.mode(M) <- "double"
  if (is.null(rownames(M))) rownames(M) <- paste0("id", seq_len(nrow(M)))

  phe <- NULL
  if (!identical(tolower(pheno_rds), "none") && file.exists(pheno_rds)) {
    phe <- readRDS(pheno_rds); if (is.matrix(phe)) phe <- as.data.frame(phe)
  }

  group <- NULL
  if (!identical(tolower(group_rds), "none") && file.exists(group_rds)) {
    gdf <- readRDS(group_rds); if (is.matrix(gdf)) gdf <- as.data.frame(gdf)
    rownames(gdf) <- as.character(gdf[[1]])
    group <- as.character(gdf[rownames(M), 2])
  }

  panel_from_matrix(M, phe, ploidy_arg, group = group, traits = traits_arg,
                    max_traits = max_traits, min_n = min_n)
}
