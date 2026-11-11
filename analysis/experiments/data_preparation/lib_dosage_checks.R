# =====================================================================
# lib_dosage_checks.R -- support checks for polyploid dosage matrices.
#
# Why this exists: the Miscanthus panel shipped for months as a MIXED-PLOIDY
# matrix analysed at a single ploidy 4.  Its source (IDB-6439486) called each
# accession at that accession's own ploidy and put everything on a common 0-2
# scale; doubling to 0-4 therefore produced genuine {0,1,2,3,4} calls for the
# tetraploid accessions but doubled-diploid {0,2,4} rows for the 2x majority.
#
# A whole-matrix support check does NOT catch that: pooled over individuals the
# support looks complete.  The individual is the unit that has a ploidy, so the
# support must be checked PER INDIVIDUAL.  On a 3-point support the single-locus
# function space is only 3-dimensional, so at ploidy 4 the degree-3 and degree-4
# orthogonal contrasts of such a row are exact linear combinations of
# {intercept, A, D2} and carry zero independent information; a panel that mixes
# those rows with real 4x rows reports kernel geometry driven by the rescaling
# rather than by autotetraploid genetics.
#
# assert_polyploid_support() is the check that was missing.  Call it in every
# build_*.R that writes a fixed-ploidy dosage matrix.
# =====================================================================

# Observed dosage classes, as a sorted integer vector.
dosage_support <- function(M) {
  v <- unique(as.vector(round(M)))
  v <- v[is.finite(v)]
  sort(as.integer(v))
}

# Rank of the single-locus function space implied by an observed support:
# k distinct dosage classes span a k-dimensional space; a ploidy-m panel needs
# m + 1 of them to make degrees 1..m independent.
support_rank <- function(support) length(support)

# Greatest common divisor of the non-zero dosage classes.  stride > 1 means the
# support is a rescaled lower-ploidy support: {0,2,4} -> 2 (doubled diploid),
# {0,3,6} -> 3 (tripled diploid), {0,1,2,3,4} -> 1 (genuine tetraploid).
support_stride <- function(support) {
  gcd2 <- function(a, b) { while (b != 0L) { t <- b; b <- a %% b; a <- t }; a }
  nz <- support[support != 0L]
  if (!length(nz)) return(1L)
  as.integer(Reduce(gcd2, nz))
}

# ---------------------------------------------------------------------
# individual_support_stride(M, max_markers)
#
# Per-row stride, i.e. one number per individual.  Subsamples columns because
# the full matrices run to hundreds of thousands of markers; with a few
# thousand markers a genuinely polyploid individual is certain to show at least
# one odd dosage class, so the subsample is not a weak test.
# ---------------------------------------------------------------------
individual_support_stride <- function(M, max_markers = 20000L, seed = 20260820L) {
  if (ncol(M) > max_markers) {
    old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
    set.seed(seed)
    M <- M[, sample.int(ncol(M), max_markers), drop = FALSE]
    if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv)
  }
  R <- round(M)
  vapply(seq_len(nrow(R)), function(i) {
    v <- R[i, ]; v <- v[is.finite(v)]
    support_stride(sort(unique(as.integer(v))))
  }, integer(1))
}

# ---------------------------------------------------------------------
# assert_polyploid_support(M, ploidy, label)
#
# Fails loudly unless the dosage matrix genuinely exercises the ploidy it
# claims.  Three hard conditions:
#   (1) values lie in 0..ploidy
#   (2) the support is NOT a rescaled lower-ploidy support, i.e. not confined
#       to multiples of some stride > 1 -- for ploidy 4 this is the {0,2,4}
#       case that motivated the check
#   (3) the support has ploidy + 1 classes, so degrees 1..ploidy are
#       identifiable rather than collinear by construction
#   (4) NO INDIVIDUAL is on a rescaled lower-ploidy support of its own -- the
#       mixed-ploidy case, which conditions (1)-(3) pass because they pool over
#       individuals
#
# min_frac guards the degenerate case where a class is technically present but
# too rare to support a contrast.
# ---------------------------------------------------------------------
assert_polyploid_support <- function(M, ploidy, label = "panel", min_frac = 1e-4,
                                     check_individuals = TRUE) {
  ploidy <- as.integer(ploidy)
  stopifnot(is.finite(ploidy), ploidy >= 2L)
  sup <- dosage_support(M)
  if (!length(sup)) stop(sprintf("[%s] dosage matrix has no finite values", label), call. = FALSE)

  if (min(sup) < 0L || max(sup) > ploidy) {
    stop(sprintf("[%s] dosage values outside 0..%d (observed %d..%d)",
                 label, ploidy, min(sup), max(sup)), call. = FALSE)
  }

  # (2) rescaled-lower-ploidy detector: is every observed class a multiple of a
  # common stride > 1?  {0,2,4} -> stride 2 -> a doubled diploid matrix.
  stride <- support_stride(sup)
  if (stride > 1L) {
    stop(sprintf(paste0("[%s] dosage support is {%s}: every class is a multiple of %d. ",
                        "This is a ploidy-%d matrix rescaled to ploidy %d, not genuine ploidy-%d ",
                        "genotype calls. Degrees above %d would be exact linear combinations of the ",
                        "lower degrees. Rebuild from calls made at ploidy %d."),
                 label, paste(sup, collapse = ","), stride,
                 as.integer(ploidy / stride), ploidy, ploidy,
                 as.integer(ploidy / stride), ploidy), call. = FALSE)
  }

  # (3) full support: every class 0..ploidy present at a usable frequency.
  want <- 0:ploidy
  frac <- vapply(want, function(k) mean(round(M) == k, na.rm = TRUE), numeric(1))
  gone <- want[frac < min_frac]
  if (length(gone)) {
    stop(sprintf(paste0("[%s] dosage classes {%s} are absent or vanishing (freq < %.1e). ",
                        "Support rank %d < %d needed for identifiable degrees 1..%d at ploidy %d."),
                 label, paste(gone, collapse = ","), min_frac,
                 support_rank(sup), ploidy + 1L, ploidy, ploidy), call. = FALSE)
  }

  # (4) per-individual support.  Conditions (1)-(3) pool over individuals and
  # so are blind to a panel that MIXES ploidies: a matrix holding real 4x rows
  # alongside doubled-diploid rows shows a complete {0..4} support overall.
  # The individual is what has a ploidy, so check the individual.
  ind_stride <- NULL
  if (check_individuals && nrow(M) > 1L) {
    ind_stride <- individual_support_stride(M)
    bad <- which(ind_stride > 1L)
    if (length(bad)) {
      ids <- rownames(M); ids <- if (is.null(ids)) as.character(bad) else ids[bad]
      stop(sprintf(paste0("[%s] %d of %d individuals are on a rescaled lower-ploidy support ",
                          "(stride > 1) while the panel is analysed at ploidy %d. ",
                          "Their implied ploidies: %s. Example IDs: %s. ",
                          "This is a MIXED-PLOIDY matrix: those rows contribute no degree-3/4 ",
                          "information and their orthogonal contrasts are exact combinations of ",
                          "the lower degrees. Rebuild from calls made at ploidy %d, or restrict ",
                          "the panel to the accessions actually called at ploidy %d."),
                   label, length(bad), nrow(M), ploidy,
                   paste(sort(unique(as.integer(ploidy / ind_stride[bad]))), collapse = "/"),
                   paste(utils::head(ids, 5), collapse = ", "), ploidy, ploidy), call. = FALSE)
    }
  }

  invisible(list(support = sup,
                 class_freq = stats::setNames(frac, want),
                 support_rank = support_rank(sup),
                 ideal_rank = ploidy + 1L,
                 individual_stride = ind_stride))
}
