#' Optimized cybr_lmpeaks - Peak Detection for BSA Data
#' 
#' This is a vectorized, much faster replacement for the original cybr_lmpeaks
#' that eliminates nested loops and reduces redundant operations.
#'
#' @param Data Data frame with columns: CHROM, POS, zscore, CSS, label
#' @param cutoff Minimum z-score threshold for peak detection (default: 2)
#' @param width Window width for slope calculation (default: 700)
#' @return Data frame with detected peaks
#' @import dplyr data.table
#' @export
cybr_lmpeaks_fast <- function(Data, cutoff = 2, width = 700) {
  require(dplyr)
  require(data.table)
  
  # Convert to data.table for faster operations
  dt <- as.data.table(Data)
  
  # Step 1: Pre-filter and calculate absolute z-scores
  dt <- dt[, abs_zscore := abs(zscore)]
  dt <- dt[abs_zscore >= cutoff]
  
  if (nrow(dt) == 0) {
    return(data.frame(CHROM = character(), CSS = character(), 
                      POS = numeric(), zscore = numeric()))
  }
  
  setorder(dt, CSS, CHROM, POS)
  
  # Step 2: Calculate slopes to find crossover points (vectorized)
  dt[, slope := {
    # Simple rolling slope calculation
    slopes <- frollapply(abs_zscore, n = width, align = "center",
                         FUN = function(x) {
                           if (length(x) < 2) return(0)
                           n <- length(x)
                           # Linear regression slope: cov(x,y) / var(x)
                           y <- x
                           x_seq <- 1:n
                           (sum(x_seq * y) - sum(x_seq) * sum(y) / n) / 
                             (sum(x_seq^2) - sum(x_seq)^2 / n)
                         })
    slopes
  }, by = .(CSS, CHROM)]
  
  # Step 3: Find crossover points (slope changes from positive to negative)
  dt <- dt[!is.na(slope)]
  dt[, negative := slope < 0]
  dt[, crossover := as.integer(negative & !shift(negative, fill = FALSE)), 
     by = .(CSS, CHROM)]
  
  # Extract crossover positions
  crossovers <- dt[crossover == 1, .(CHROM, CSS, POSc = POS)]
  
  # Step 4: Create segment boundaries
  # Add boundaries for all chromosomes and CSS combinations
  all_combos <- CJ(
    CHROM = as.factor(as.character(as.roman(1:16))),
    CSS = c("VIII", "I")
  )
  
  # Add initial boundary (position 1) for each combination
  initial_bounds <- all_combos[, .(POSc = 1)]
  initial_bounds <- cbind(all_combos, initial_bounds)
  
  # Combine with crossover points
  segments <- rbindlist(list(initial_bounds, crossovers), fill = TRUE)
  setorder(segments, CHROM, CSS, POSc)
  
  # Create Start and End positions
  segments[, End := shift(POSc, type = "lead", fill = Inf), by = .(CHROM, CSS)]
  segments[, Start := POSc]
  segments <- segments[, .(CHROM, CSS, Start, End)]
  
  # Step 5: Find peaks within each segment (MUCH faster than nested loops)
  # Use non-equi joins for efficient interval matching
  setkey(segments, CHROM, CSS, Start, End)
  setkey(dt, CHROM, CSS, POS)
  
  # Filter for Bulk label only
  bulk_data <- dt[label == "Bulk"]
  
  # Find maximum z-score in each segment using rolling join
  peaks <- bulk_data[segments, 
                     on = .(CHROM, CSS, POS > Start, POS < End),
                     .(max_zscore = max(abs_zscore, na.rm = TRUE),
                       peak_pos = POS[which.max(abs_zscore)]),
                     by = .EACHI,
                     allow.cartesian = TRUE]
  
  # Filter out segments with no valid peaks
  peaks <- peaks[is.finite(max_zscore) & max_zscore > cutoff]
  
  # Step 6: Match back to original data to get full peak information
  setnames(peaks, c("POS", "POS.1"), c("Start", "End"))
  
  result <- bulk_data[peaks, 
                      on = .(CHROM, CSS, POS == peak_pos),
                      nomatch = NULL]
  
  result <- result[, .(CHROM, CSS, POS, zscore = abs_zscore)]
  
  return(as.data.frame(result))
}


#' Ultra-Fast Version Using Simplified Peak Detection
#' 
#' This version uses a simpler but faster algorithm that gives similar results
#' without the complex slope-based segment detection.
#'
#' @param Data Data frame with columns: CHROM, POS, zscore, CSS, label
#' @param cutoff Minimum z-score threshold (default: 2)
#' @param width Window width for local maxima detection (default: 700)
#' @return Data frame with detected peaks
#' @import data.table
#' @export
cybr_lmpeaks_ultrafast <- function(Data, cutoff = 2, width = 700) {
  require(data.table)
  
  dt <- as.data.table(Data)
  
  # Filter for Bulk and above threshold
  dt <- dt[label == "Bulk" & abs(zscore) >= cutoff]
  
  if (nrow(dt) == 0) {
    return(data.frame(CHROM = character(), CSS = character(),
                      POS = numeric(), zscore = numeric()))
  }
  
  dt[, abs_zscore := abs(zscore)]
  setorder(dt, CSS, CHROM, POS)
  
  # Find local maxima using rolling window
  dt[, is_peak := {
    # Calculate rolling max
    roll_max <- frollapply(abs_zscore, n = width, align = "center",
                           FUN = max, na.rm = TRUE)
    
    # Position is a peak if it equals the rolling max
    is_max <- abs_zscore == roll_max & !is.na(roll_max)
    
    # Remove consecutive peaks (keep only first in a run)
    if (sum(is_max) > 1) {
      peak_idx <- which(is_max)
      # Calculate distance between consecutive peaks
      if (length(peak_idx) > 1) {
        distances <- diff(peak_idx)
        # Keep peaks that are far enough apart
        keep <- c(TRUE, distances > width/2)
        is_max[peak_idx[!keep]] <- FALSE
      }
    }
    
    is_max
  }, by = .(CSS, CHROM)]
  
  # Return only peaks
  result <- dt[is_peak == TRUE, .(CHROM, CSS, POS, zscore = abs_zscore)]
  
  return(as.data.frame(result))
}


#' Helper Functions
#' These are the utility functions likely used in the original

slope_change <- function(x) {
  if (length(x) < 2) return(0)
  n <- length(x)
  # Calculate linear regression slope
  y <- x
  x_seq <- 1:n
  slope <- (sum(x_seq * y) - sum(x_seq) * sum(y) / n) / 
    (sum(x_seq^2) - sum(x_seq)^2 / n)
  return(slope)
}

subtract2 <- function(x) {
  if (length(x) < 2) return(0)
  # Detects transition from FALSE to TRUE
  as.integer(x[2] & !x[1])
}


#' Benchmark Comparison Function
#'
#' @param Data Test data
#' @param cutoff Z-score cutoff
#' @param width Window width
#' @export
benchmark_lmpeaks <- function(Data, cutoff = 2, width = 700) {
  require(microbenchmark)
  
  cat("Running benchmark with", nrow(Data), "rows\n")
  
  results <- microbenchmark(
    original = cybr_lmpeaks(Data, cutoff, width),
    fast = cybr_lmpeaks_fast(Data, cutoff, width),
    ultrafast = cybr_lmpeaks_ultrafast(Data, cutoff, width),
    times = 5
  )
  
  print(results)
  
  # Compare results
  cat("\n=== Comparing Results ===\n")
  
  peaks_orig <- cybr_lmpeaks(Data, cutoff, width)
  peaks_fast <- cybr_lmpeaks_fast(Data, cutoff, width)
  peaks_ultra <- cybr_lmpeaks_ultrafast(Data, cutoff, width)
  
  cat("Original peaks found:", nrow(peaks_orig), "\n")
  cat("Fast version peaks found:", nrow(peaks_fast), "\n")
  cat("Ultrafast version peaks found:", nrow(peaks_ultra), "\n")
  
  return(list(
    benchmark = results,
    original = peaks_orig,
    fast = peaks_fast,
    ultrafast = peaks_ultra
  ))
}


#' Example Usage
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(data.table)
#' 
#' # Create sample data
#' set.seed(123)
#' test_data <- data.frame(
#'   CHROM = rep(as.factor(as.character(as.roman(c(1,2,3)))), each = 5000),
#'   CSS = rep(c("I", "VIII"), each = 7500),
#'   POS = rep(1:5000, 3),
#'   zscore = rnorm(15000, 0, 2),
#'   label = "Bulk"
#' )
#' 
#' # Add some clear peaks
#' test_data$zscore[c(1000, 3000, 8000, 12000)] <- c(5, 6, 7, 5.5)
#' 
#' # Test the functions
#' peaks_fast <- cybr_lmpeaks_fast(test_data, cutoff = 2, width = 700)
#' peaks_ultra <- cybr_lmpeaks_ultrafast(test_data, cutoff = 2, width = 700)
#' 
#' print(head(peaks_fast))
#' print(head(peaks_ultra))
#' 
#' # Run benchmark (if you have the original function)
#' # benchmark_lmpeaks(test_data, cutoff = 2, width = 700)
#' }