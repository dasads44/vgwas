#' Perform Association Studies Across Genomic Data
#'
#' @description
#' This function performs association studies across genomic data, supporting two different modes:
#' 1) comparing each amino acid with count > aa.freq.cutoff against all others at each site, or
#' 2) comparing only the most common amino acid against all others at each site.
#'
#' @param data A data.frame or data.table containing rows per sample with required model data
#'        and amino acids for each position in the genome labeled X1:Xn.
#' @param start Integer. Amino acid position to start the analysis (default = 1).
#' @param end Integer. Amino acid position to end the analysis.
#' @param formula Formula object. Model formula to use for logistic regression on each site.
#' @param aa.freq.cutoff Integer. Minimum number of samples required for an amino acid to be included in test.
#' @param fdr Numeric. False discovery rate threshold between 0 and 1.
#' @param mode Character. Either "all" (compare each amino acid against all others) or
#'        "dominant" (compare only the most common amino acid against all others).
#'
#' @return A data.table containing the association study results ordered by p-value.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data <- read.csv("genome_data.csv")
#'
#' # Define formula for model
#' model_formula <- outcome ~ site + covariate1 + covariate2
#'
#' # Run analysis comparing all amino acids at each position
#' results_all <- genomeAssociation(
#'   data = data,
#'   start = 1,
#'   end = 500,
#'   formula = model_formula,
#'   aa.freq.cutoff = 10,
#'   fdr = 0.05,
#'   mode = "all"
#' )
#'
#' # Run analysis comparing only dominant amino acid
#' results_dominant <- genomeAssociation(
#'   data = data,
#'   start = 1,
#'   end = 500,
#'   formula = model_formula,
#'   aa.freq.cutoff = 10,
#'   fdr = 0.05,
#'   mode = "dominant"
#' )
#' }
#'
#' @importFrom qvalue qvalue
#' @importFrom data.table as.data.table
#' @export
genomeAssociation <- function(data,
                              start = 1,
                              end,
                              formula,
                              aa.freq.cutoff,
                              fdr,
                              mode = c("all", "dominant")) {
  # Check inputs
  mode <- match.arg(mode)

  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame or data.table")
  }

  if (!is.numeric(start) || !is.numeric(end) || start < 1 || end < start) {
    stop("Invalid 'start' or 'end' parameters")
  }

  if (!is.numeric(aa.freq.cutoff) || aa.freq.cutoff < 1) {
    stop("'aa.freq.cutoff' must be a positive integer")
  }

  if (!is.numeric(fdr) || fdr <= 0 || fdr >= 1) {
    stop("'fdr' must be between 0 and 1")
  }

  # Converting input data to data.frame
  data <- as.data.frame(data)
  temp.df <- as.data.frame(data)

  # Initialize results data frame
  results <- data.frame(site = integer(),
                        aa = character(),
                        p.val = double(),
                        log_OR = double(),
                        SE = double())

  # Iterate through genomic positions
  for (temp.Var in paste('X', start:end, sep = '')) {
    # Skip if column doesn't exist
    if (!(temp.Var %in% names(data))) {
      next
    }

    temp.idx <- which(names(data) == temp.Var)
    temp.aa <- sort(table(data[, temp.idx]), decreasing = TRUE)

    if (mode == "all") {
      # Mode 1: Compare each amino acid with sufficient count against all others
      if (length(temp.aa) == 2) {
        if (temp.aa[2] >= aa.freq.cutoff) {
          analyzePosition(temp.Var, temp.idx, temp.aa, 1, data, temp.df, formula, results)
        }
      } else if (length(temp.aa) >= 3) {
        for (i in 1:length(temp.aa)) {
          if (temp.aa[i] >= aa.freq.cutoff & sum(temp.aa[-i]) >= aa.freq.cutoff) {
            results <- analyzePosition(temp.Var, temp.idx, temp.aa, i, data, temp.df, formula, results)
          }
        }
      }
    } else if (mode == "dominant") {
      # Mode 2: Compare only the most common amino acid against all others
      if (length(temp.aa) >= 2) {
        if (sum(temp.aa[-1]) >= aa.freq.cutoff) {
          results <- analyzePosition(temp.Var, temp.idx, temp.aa, 1, data, temp.df, formula, results)
        }
      }
    }
  }

  # Early return if no results
  if (nrow(results) == 0) {
    warning("No significant associations found with the given parameters")
    return(data.table::data.table())
  }

  # Calculating FDR values
  fdr.res <- qvalue::qvalue(results$p.val, fdr.level = fdr)
  fdr.value <- -log10(max(fdr.res$pvalues[fdr.res$significant == TRUE]))

  # Convert to data.table for better performance
  results <- data.table::as.data.table(results)

  # Add q-values and FDR threshold
  results$q.val <- fdr.res$qvalues

  # Return results ordered by p-value
  return(list(results = results[order(p.val)],
              data = data,
              log10.fdr.value = fdr.value))
}

#' Helper function to analyze a single position
#'
#' @param temp.Var Current variable name
#' @param temp.idx Index of current variable in data
#' @param temp.aa Table of amino acid frequencies
#' @param aa_index Index of the amino acid being analyzed
#' @param data Original data
#' @param temp.df Working data frame
#' @param formula Formula for the model
#' @param results Results data frame to append to
#'
#' @return Updated results data frame
#' @keywords internal
analyzePosition <- function(temp.Var, temp.idx, temp.aa, aa_index, data, temp.df, formula, results) {
  # Create factor for the current amino acid vs all others
  temp.df$site <- factor(x = data[, temp.idx] == names(temp.aa)[aa_index],
                         levels = c(FALSE, TRUE),
                         labels = c(paste0('non_', names(temp.aa)[aa_index]),
                                    names(temp.aa)[aa_index]))

  temp.df$site <- stats::relevel(x = temp.df$site,
                                 ref = paste0('non_', names(temp.aa)[aa_index]))

  # Fit logistic regression model
  temp.model <- stats::glm(formula = formula,
                           data = temp.df,
                           family = stats::binomial)

  # Extract coefficients
  tt <- stats::coef(summary(temp.model))

  # Check if the model fit correctly
  if (nrow(tt) < 2 || ncol(tt) < 4) {
    return(results)  # Return unchanged if model didn't fit correctly
  }

  # Create result row
  temp.res <- data.frame(site = as.numeric(gsub('X', '', temp.Var)),
                         aa = names(temp.aa)[aa_index],
                         p.val = tt[2, 4],
                         log_OR = tt[2, 1],
                         SE = tt[2, 2])

  # Append to results
  return(rbind(results, temp.res))
}
