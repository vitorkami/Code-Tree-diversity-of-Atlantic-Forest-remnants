# =============================================================================
# fitoR: Phytosociological descriptors
#
# Adapted/translated to English for reproducible workflows.
#
# Original repository: https://github.com/ricds/fitoR
# Author: Ricardo Dal'Agnol
# Contributors (as listed in the original header): Alexandre Gabriel Christo,
# Pedro Higuchi, Arthur Vinicius Rodrigues
#
# Notes (from original documentation, translated):
# - Input data must include:
#   * parc: plot/parcel identifier (character or numeric)
#   * spp : species/taxon identifier
#   * dap (diameter at breast height) OR cap (circumference at breast height)
# - If an individual has multiple stems, create dap1, dap2, ... (or cap1, cap2, ...)
# - The function searches for columns containing 'dap' or 'cap'.
#   Do not include unrelated columns containing these strings (e.g., 'dap_mean'),
#   otherwise basal area may be miscomputed.
#
# Function:
#   fitoR(x, area_m2_per_plot, filename = NULL)
# =============================================================================

fitoR <- function(x, area_m2_per_plot, filename = NULL) {
  
  stopifnot(all(c("parc", "spp") %in% colnames(x)))
  x <- as.data.frame(x)
  
  # Replace NA with 0 only in numeric columns (safer than x[is.na(x)] <- 0)
  num_cols <- vapply(x, is.numeric, logical(1))
  x[num_cols][is.na(x[num_cols])] <- 0
  
  # Species-by-plot abundance table (counts of individuals)
  mat <- table(x$spp, x$parc)
  
  n_plots <- length(unique(x$parc))
  sampled_area_m2 <- area_m2_per_plot * n_plots
  sampled_area_ha <- sampled_area_m2 / 10000
  
  # Individuals per species
  N <- apply(mat, 1, sum)
  
  # Density (absolute and relative)
  DA <- N / sampled_area_ha
  DR <- (DA / sum(DA)) * 100
  
  # Frequency (absolute and relative)
  freq <- if (length(dim(mat)) > 1) apply(mat > 0, 1, sum) else sum(mat > 0)
  FA <- (freq / n_plots) * 100
  FR <- (FA / sum(FA)) * 100
  
  # Detect diameter/circumference columns
  cap_cols <- grep("cap", colnames(x), ignore.case = TRUE)
  dap_cols <- grep("dap", colnames(x), ignore.case = TRUE)
  
  if (length(dap_cols) > 0) {
    param <- "dap"
    cols <- dap_cols
  } else if (length(cap_cols) > 0) {
    param <- "cap"
    cols <- cap_cols
  } else {
    stop("No 'dap' or 'cap' columns found in input data.")
  }
  
  # Basal area per individual (m2) summed across multiple stems
  # If CAP: basal area = pi*(d^2)/4 with d = cap/pi
  # If DAP: basal area = pi*(dap^2)/4
  # Convert cm^2 to m^2: divide by 10000
  x$basal_area_m2 <- 0
  
  for (i in cols) {
    if (param == "cap") {
      d <- x[[i]] / pi
      x$basal_area_m2 <- x$basal_area_m2 + (pi * (d^2) / 4) / 10000
    }
    if (param == "dap") {
      d <- x[[i]]
      x$basal_area_m2 <- x$basal_area_m2 + (pi * (d^2) / 4) / 10000
    }
  }
  
  # Dominance (absolute and relative)
  DoA <- tapply(x$basal_area_m2, x$spp, sum) / sampled_area_ha
  DoR <- (DoA / sum(DoA)) * 100
  
  # Total basal area per species (m2)
  BA <- tapply(x$basal_area_m2, x$spp, sum)
  
  # Importance value index (mean of DR, FR, DoR)
  VI <- (DR + FR + DoR) / 3
  
  fito <- data.frame(
    N = as.numeric(N),
    BA = as.numeric(BA),
    DA = as.numeric(DA),
    DR = as.numeric(DR),
    DoA = as.numeric(DoA),
    DoR = as.numeric(DoR),
    FA = as.numeric(FA),
    FR = as.numeric(FR),
    VI = as.numeric(VI)
  )
  
  # Rounding + ordering
  fito$BA  <- round(fito$BA, 2)
  fito$DA  <- round(fito$DA, 2)
  fito$DR  <- round(fito$DR, 2)
  fito$DoA <- round(fito$DoA, 2)
  fito$DoR <- round(fito$DoR, 2)
  fito$FA  <- round(fito$FA, 2)
  fito$FR  <- round(fito$FR, 2)
  fito$VI  <- round(fito$VI, 2)
  
  fito <- fito[order(fito$VI, decreasing = TRUE), , drop = FALSE]
  
  if (!is.null(filename)) {
    utils::write.csv(fito, file = paste0(filename, ".csv"))
  }
  
  return(fito)
}
