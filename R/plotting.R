#' Generate a Quantile-Quantile Plot for Association Results
#'
#' @description
#' Creates a quantile-quantile (QQ) plot to visualize the distribution of p-values
#' from genome-wide association studies. This helps to assess whether the observed
#' p-values deviate from the expected distribution under the null hypothesis.
#'
#' @param results A data.frame or data.table output from the genomeAssociation function,
#'                containing at minimum a 'p.val' column with association p-values.
#'
#' @details
#' The QQ plot compares the observed -log10(p-values) against their expected values
#' under the null hypothesis (uniform distribution). Deviations from the diagonal line
#' indicate potential systematic bias in the test statistics or the presence of true
#' associations.
#'
#' The x-axis represents the expected -log10(p) values, calculated as:
#' -log10(rank/(n+1)), where rank is the rank of each p-value and n is the total number
#' of tests.
#'
#' The y-axis represents the observed -log10(p) values.
#'
#' @return A ggplot2 object containing the QQ plot, which can be further modified or
#'         directly printed.
#'
#' @examples
#' \dontrun{
#' # Generate QQ plot from association results
#' data <- read.csv("genome_data.csv")
#' formula <- outcome ~ site + covariate
#' results <- genomeAssociation(data, 1, 500, formula, 10, 0.05, mode="all")
#' qq_plot <- qqggplot(results)
#'
#' # Display the plot
#' print(qq_plot)
#'
#' # Modify the plot if needed
#' qq_plot +
#'   theme_minimal() +
#'   ggtitle("QQ Plot of Genomic Association Results")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline labs
#' @export
qqggplot <- function(results){
  # Check input
  if (!("p.val" %in% names(results))) {
    stop("Input 'results' must contain a 'p.val' column")
  }

  if (!is.numeric(results$p.val)) {
    stop("The 'p.val' column must contain numeric values")
  }

  if (any(results$p.val < 0 | results$p.val > 1, na.rm = TRUE)) {
    warning("Some p-values are outside the range [0,1]. Check your input data.")
  }

  # Remove NA values if present
  results <- results[!is.na(results$p.val), ]

  # Create QQ plot
  qqplot <- ggplot2::ggplot(data = results,
                            ggplot2::aes(x = -log10(rank(p.val, ties.method = "first")/(length(p.val) + 1)),
                                         y = -log10(p.val))) +
    ggplot2::geom_point(shape = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
    ggplot2::labs(y = "Observed -log10(p)",
                  x = "Expected -log10(p)",
                  title = "Quantile-Quantile Plot of Association p-values")

  return(qqplot)
}

#' Plot Manhattan-style GWAS Results
#'
#' @description
#' Creates a Manhattan-style plot for genome-wide association study (GWAS) results,
#' highlighting significant sites that exceed the FDR threshold. The function
#' selects the most significant p-value for each site and labels significant points.
#'
#' @param results A list containing two elements:
#'   \itemize{
#'     \item \code{results}: A data.table or data.frame with GWAS results. Must contain columns:
#'       \itemize{
#'         \item \code{site}: Genomic position or site identifier
#'         \item \code{p.val}: P-values for each site
#'         \item \code{aa}: Amino acid or variant information (optional)
#'       }
#'     \item \code{log10.fdr.value}: Numeric value representing -log10 of the FDR threshold
#'   }
#' @param max.num.overlaps Integer specifying the maximum number of overlapping labels to display.
#'   Higher values show more labels but may cause crowding. Default is 20.
#'
#' @return A ggplot2 object showing the Manhattan plot of GWAS results
#'
#' @importFrom ggplot2 ggplot aes labs geom_point geom_hline scale_y_continuous scale_x_continuous
#'   theme_bw theme element_blank element_text scale_color_manual
#' @importFrom ggrepel geom_label_repel
#' @importFrom data.table .SD is.data.table as.data.table
#'
#' @examples
#' \dontrun{
#' # Example data structure
#' gwas_data <- list(
#'   results = data.table(
#'     site = rep(1:100, each = 3),
#'     p.val = runif(300, 0, 0.1),
#'     aa = sample(LETTERS, 300, replace = TRUE)
#'   ),
#'   log10.fdr.value = 2
#' )
#'
#' # Create and display the plot
#' Plot_GWAS_Fig(gwas_data)
#'
#' # Adjust maximum number of overlapping labels
#' Plot_GWAS_Fig(gwas_data, max.num.overlaps = 10)
#' }
#'
Plot_GWAS_Fig <- function(results, max.num.overlaps = 20) {
  # Input validation
  if (!is.list(results) || !all(c("results", "log10.fdr.value") %in% names(results))) {
    stop("Input must be a list with 'results' and 'log10.fdr.value' elements")
  }

  # Extract components from the results list
  fdr.val <- results$log10.fdr.value
  result_data <- results$results

  # Check if result_data is a data.frame or data.table
  if (!is.data.frame(result_data)) {
    stop("The 'results' element must be a data.frame or data.table")
  }

  # Required columns check
  required_cols <- c("site", "p.val")
  if (!all(required_cols %in% colnames(result_data))) {
    stop(paste("Missing required columns in results:",
               paste(required_cols[!required_cols %in% colnames(result_data)], collapse = ", ")))
  }

  # Convert to data.table if it's not already
  if (!requireNamespace("data.table", quietly = TRUE)) {
    if (!is.data.table(result_data)) {
      stop("Package 'data.table' is required but not installed")
    }
  } else if (!data.table::is.data.table(result_data)) {
    result_data <- data.table::as.data.table(result_data)
  }

  # Check if FDR value is valid
  if (is.nan(fdr.val)) {
    warning("No significant results found. FDR cutoff line will not be displayed.")
  }

  # Get the most significant p-value for each site
  plot_data <- result_data[, .SD[which.min(p.val)], by = site]

  # Create the plot
  plot <- ggplot(data = plot_data,
                 aes(x = as.numeric(site), y = -log10(p.val))) +
    labs(y = expression(-log[10]~p), x = "Genomic Position") +
    geom_point(aes(color = -log10(p.val) >= fdr.val), size = 1.2) +
    geom_hline(yintercept = fdr.val, linetype = "dashed", color = "darkgray") +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(0, max(as.numeric(plot_data$site))),
                       expand = c(0.01, 0)) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red3")) +
    geom_label_repel(
      aes(label = ifelse(test = -log10(p.val) >= fdr.val,
                         yes = paste0(site, toupper(aa)),
                         no = '')),
      colour = "black",
      point.padding = unit(0.5, "lines"),
      box.padding = unit(0.5, 'lines'),
      arrow = arrow(length = unit(0.01, 'npc')),
      na.rm = TRUE,
      max.overlaps = max.num.overlaps
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold")
    )

  return(plot)
}
