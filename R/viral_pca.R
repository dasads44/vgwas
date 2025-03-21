#' Convert Sequence Alignment to Binary Matrix for PCA
#'
#' This function converts a sequence alignment object (as returned by \code{seqinr::read.alignment})
#' into a binary matrix suitable for principal component analysis (PCA). Each column in the
#' resulting matrix represents the presence (1) or absence (0) of a specific amino acid or
#' nucleotide at a specific position in the alignment.
#'
#' @param alignment An alignment object with the following components:
#'   \describe{
#'     \item{seq}{A list of character strings representing sequences}
#'     \item{nam}{A character vector of sequence names}
#'   }
#' @param exclude_gaps Logical. If TRUE, gaps ("-") are excluded from the binary encoding. Default is TRUE.
#' @param min_freq Numeric between 0 and 1. Minimum frequency for a variant to be included in the binary matrix.
#'   Variants with frequency below this threshold will be excluded. Default is 0 (include all variants).
#' @param verbose Logical. If TRUE, progress messages are printed. Default is FALSE.
#'
#' @return A binary matrix where rows correspond to sequences and columns represent the presence
#'   or absence of specific amino acids or nucleotides at specific positions. Column names are
#'   formatted as "posX_Y" where X is the position and Y is the amino acid or nucleotide.
#'
#' @examples
#' \dontrun{
#' # Read an alignment file
#' library(seqinr)
#' my_alignment <- read.alignment(file = "my_alignment.fasta", format = "fasta")
#'
#' # Convert to binary matrix
#' bin_matrix <- alignment_to_binary_matrix(my_alignment)
#'
#' # Perform PCA
#' pca_result <- prcomp(bin_matrix, scale = TRUE)
#'
#' # Plot the first two principal components
#' plot(pca_result$x[,1], pca_result$x[,2],
#'      xlab = "PC1", ylab = "PC2",
#'      main = "PCA of sequence alignment")
#' }
#'
#' @export
#' @importFrom stats prcomp
alignment_to_binary_matrix <- function(alignment, exclude_gaps = TRUE, min_freq = 0, verbose = FALSE) {
  # Input validation
  if (!is.list(alignment) || is.null(alignment$seq) || is.null(alignment$nam)) {
    stop("Alignment must be a list with 'seq' and 'nam' components")
  }

  if (min_freq < 0 || min_freq > 1) {
    stop("min_freq must be between 0 and 1")
  }

  # Get the number of sequences and sites
  n_seq <- length(alignment$seq)
  seq_length <- nchar(alignment$seq[[1]])

  if (verbose) {
    message(sprintf("Processing alignment with %d sequences of length %d", n_seq, seq_length))
  }

  # Use lapply to split all sequences at once
  if (verbose) message("Splitting sequences into characters...")
  seq_chars <- lapply(alignment$seq, function(seq) strsplit(seq, "")[[1]])

  # Fill the character matrix
  char_matrix <- do.call(rbind, seq_chars)

  if (verbose) message("Creating binary matrix...")

  # Process each position to create binary columns
  position_cols <- lapply(1:seq_length, function(j) {
    if (verbose && j %% 100 == 0) {
      message(sprintf("Processing position %d of %d", j, seq_length))
    }

    # Get unique amino acids at this position
    unique_aa <- unique(char_matrix[, j])

    # Exclude gaps if requested
    if (exclude_gaps) {
      unique_aa <- unique_aa[unique_aa != "-"]
    }

    # Filter by minimum frequency if specified
    if (min_freq > 0) {
      aa_freq <- sapply(unique_aa, function(aa) {
        sum(char_matrix[, j] == aa) / n_seq
      })
      unique_aa <- unique_aa[aa_freq >= min_freq]
    }

    # If no variants remain after filtering, return NULL
    if (length(unique_aa) == 0) {
      return(NULL)
    }

    # Create binary matrix for this position
    pos_matrix <- sapply(unique_aa, function(aa) {
      as.integer(char_matrix[, j] == aa)
    })

    # Set column names
    colnames(pos_matrix) <- paste0("pos", j, "_", unique_aa)

    return(pos_matrix)
  })

  # Remove NULL elements (positions with no variants passing filters)
  position_cols <- position_cols[!sapply(position_cols, is.null)]

  # Combine all position matrices into one
  if (length(position_cols) == 0) {
    stop("No variants passed the filtering criteria. Try lowering min_freq.")
  }

  binary_matrix <- do.call(cbind, position_cols)
  rownames(binary_matrix) <- alignment$nam

  if (verbose) {
    message(sprintf("Binary matrix created with %d rows and %d columns",
                    nrow(binary_matrix), ncol(binary_matrix)))
  }

  return(binary_matrix)
}

#' Perform PCA on Sequence Alignment
#'
#' This function converts a sequence alignment to a binary matrix and performs
#' principal component analysis (PCA) on it.
#'
#' @param alignment An alignment object with 'seq' and 'nam' components
#' @param exclude_gaps Logical. If TRUE, gaps ("-") are excluded from the binary encoding. Default is TRUE.
#' @param min_freq Numeric between 0 and 1. Minimum frequency for a variant to be included. Default is 0.
#' @param center Logical. If TRUE, the variables are centered before PCA. Default is TRUE.
#' @param scale Logical. If TRUE, the variables are scaled before PCA. Default is TRUE.
#' @param verbose Logical. If TRUE, progress messages are printed. Default is FALSE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{pca}{The result of prcomp}
#'     \item{binary_matrix}{The binary matrix used for PCA}
#'     \item{var_explained}{The proportion of variance explained by each principal component}
#'   }
#'
#' @examples
#' \dontrun{
#' # Read an alignment file
#' library(seqinr)
#' my_alignment <- read.alignment(file = "my_alignment.fasta", format = "fasta")
#'
#' # Perform PCA
#' pca_result <- alignment_pca(my_alignment, min_freq = 0.01)
#'
#' # Plot the first two principal components
#' plot(pca_result$pca$x[,1], pca_result$pca$x[,2],
#'      xlab = "PC1", ylab = "PC2",
#'      main = "PCA of sequence alignment")
#'
#' # Plot variance explained
#' barplot(pca_result$var_explained[1:10],
#'         names.arg = paste0("PC", 1:10),
#'         main = "Scree Plot",
#'         ylab = "Proportion of Variance Explained")
#' }
#'
#' @export
#' @importFrom stats prcomp
alignment_pca <- function(alignment, exclude_gaps = TRUE, min_freq = 0,
                          center = TRUE, scale = TRUE, verbose = FALSE) {
  # Convert alignment to binary matrix
  binary_matrix <- alignment_to_binary_matrix(alignment,
                                              exclude_gaps = exclude_gaps,
                                              min_freq = min_freq,
                                              verbose = verbose)

  if (verbose) message("Performing PCA...")

  # Perform PCA
  pca_result <- prcomp(binary_matrix, center = center, scale = scale)

  # Calculate variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

  # Return results
  return(list(
    pca = pca_result,
    binary_matrix = binary_matrix,
    var_explained = var_explained
  ))
}

#' Plot PCA Results from Sequence Alignment
#'
#' This function creates a scatter plot of the first two principal components
#' from a PCA analysis of sequence alignment data.
#'
#' @param pca_result The result from \code{alignment_pca} or \code{prcomp}
#' @param groups An optional factor or vector for coloring points by group
#' @param pc_x Integer. The principal component to plot on the x-axis. Default is 1.
#' @param pc_y Integer. The principal component to plot on the y-axis. Default is 2.
#' @param main Character. The main title for the plot. Default is "PCA of sequence alignment".
#' @param ... Additional arguments passed to \code{plot}
#'
#' @return Invisibly returns the plot coordinates
#'
#' @examples
#' \dontrun{
#' # Read an alignment file
#' library(seqinr)
#' my_alignment <- read.alignment(file = "my_alignment.fasta", format = "fasta")
#'
#' # Perform PCA
#' pca_result <- alignment_pca(my_alignment)
#'
#' # Basic plot
#' plot_alignment_pca(pca_result)
#'
#' # Plot with grouping
#' genotypes <- factor(c("A", "B", "A", ...))  # One label per sequence
#' plot_alignment_pca(pca_result, groups = genotypes,
#'                   main = "PCA by Genotype")
#' }
#'
#' @export
#' @importFrom graphics plot legend
plot_alignment_pca <- function(pca_result, groups = NULL, pc_x = 1, pc_y = 2,
                               main = "PCA of sequence alignment", ...) {
  # Extract PCA object if nested in alignment_pca result
  if (is.list(pca_result) && !is.null(pca_result$pca)) {
    pca_obj <- pca_result$pca
    var_explained <- pca_result$var_explained
  } else {
    pca_obj <- pca_result
    var_explained <- pca_obj$sdev^2 / sum(pca_obj$sdev^2)
  }

  # Check if principal components exist
  if (pc_x > ncol(pca_obj$x) || pc_y > ncol(pca_obj$x)) {
    stop("Requested principal components exceed the number available")
  }

  # Create axis labels with variance explained
  xlab <- sprintf("PC%d (%.1f%%)", pc_x, var_explained[pc_x] * 100)
  ylab <- sprintf("PC%d (%.1f%%)", pc_y, var_explained[pc_y] * 100)

  # Prepare plot parameters
  plot_params <- list(
    x = pca_obj$x[, pc_x],
    y = pca_obj$x[, pc_y],
    xlab = xlab,
    ylab = ylab,
    main = main,
    pch = 19,
    ...
  )

  # Add color by group if provided
  if (!is.null(groups)) {
    if (length(groups) != nrow(pca_obj$x)) {
      stop("Length of groups must match the number of sequences")
    }

    groups <- factor(groups)
    colors <- rainbow(length(levels(groups)))
    plot_params$col <- colors[groups]

    # Create the plot
    do.call(plot, plot_params)

    # Add legend
    legend("topright", legend = levels(groups), col = colors, pch = 19)
  } else {
    # Create the plot without grouping
    do.call(plot, plot_params)
  }

  # Return coordinates invisibly
  invisible(list(x = pca_obj$x[, pc_x], y = pca_obj$x[, pc_y]))
}
