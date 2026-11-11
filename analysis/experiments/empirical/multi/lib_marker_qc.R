marker_qc_enabled <- function() Sys.getenv("PWF_MARKER_QC", "1") != "0"

marker_qc_maf_min <- function() {
  x <- suppressWarnings(as.numeric(Sys.getenv("PWF_MAF_MIN", "0.05")))
  if (!is.finite(x) || x < 0) stop("PWF_MAF_MIN must be a non-negative number", call. = FALSE)
  x
}

marker_qc_call_min <- function() {
  x <- suppressWarnings(as.numeric(Sys.getenv("PWF_CALLRATE_MIN", "0.8")))
  if (!is.finite(x) || x < 0 || x > 1) stop("PWF_CALLRATE_MIN must be in [0, 1]", call. = FALSE)
  x
}

apply_marker_qc <- function(M, ploidy, marker_ploidy = NULL,
                            enabled = marker_qc_enabled(),
                            maf_min = marker_qc_maf_min(),
                            call_min = marker_qc_call_min()) {
  if (is.data.frame(M)) M <- as.matrix(M)
  storage.mode(M) <- "double"
  n0 <- nrow(M)
  m0 <- ncol(M)
  if (!is.null(marker_ploidy)) {
    marker_ploidy <- as.integer(marker_ploidy)
    if (length(marker_ploidy) != m0) {
      stop("marker ploidy vector length must equal ncol(M)", call. = FALSE)
    }
    denom <- marker_ploidy
    ploidy_label <- "variable"
  } else {
    ploidy <- as.integer(ploidy)
    if (!is.finite(ploidy) || ploidy < 2) stop("fixed ploidy must be >=2", call. = FALSE)
    denom <- rep(ploidy, m0)
    ploidy_label <- as.character(ploidy)
  }

  missing_pre_qc <- sum(is.na(M))
  callrate <- colMeans(!is.na(M))
  af <- colMeans(M, na.rm = TRUE) / denom
  maf <- pmin(af, 1 - af)
  polymorphic <- apply(M, 2, function(z) length(unique(na.omit(round(z)))) > 1)
  keep <- is.finite(callrate) & callrate > call_min & is.finite(maf) &
    maf >= maf_min & polymorphic
  # Keep every *usable* marker; all-missing markers can never be imputed, and
  # letting them through leaves NaN cells that only fail later. Must stay in
  # step with polyDx::filter_markers(), which writes the same markers_post_qc.
  if (!enabled) keep <- is.finite(callrate) & callrate > 0
  if (!any(keep)) stop("no markers remain after marker QC", call. = FALSE)

  M <- M[, keep, drop = FALSE]
  if (!is.null(marker_ploidy)) marker_ploidy <- marker_ploidy[keep]
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
  M[M < 0] <- 0
  if (is.null(marker_ploidy)) {
    M[M > ploidy] <- ploidy
  } else {
    MP <- matrix(marker_ploidy, nrow(M), ncol(M), byrow = TRUE)
    too_high <- M > MP
    M[too_high] <- MP[too_high]
  }
  storage.mode(M) <- "double"

  qc <- data.frame(
    ploidy = ploidy_label,
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
    median_maf_pre_qc = suppressWarnings(median(maf, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  list(M = M, marker_ploidy = marker_ploidy, qc = qc)
}

write_marker_qc_audit <- function(prefix, dataset, qc) {
  out <- cbind(data.frame(dataset = dataset, stringsAsFactors = FALSE), qc)
  tmp <- paste0(prefix, "_qcaudit.csv.tmp", Sys.getpid())
  write.csv(out, tmp, row.names = FALSE)
  file.rename(tmp, paste0(prefix, "_qcaudit.csv"))
}
