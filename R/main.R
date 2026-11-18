# Functions reused from Practical 5. Add documentation in your new package.
#' Encode DNA sequences as position-specific factors
#'
#' Split DNA sequences into individual nucleotides and represent each sequence
#' position as a categorical factor for use by the prediction model. This is an
#' internal helper; it does not perform one-hot encoding.
#'
#' @param dna_strings A non-empty character vector of non-missing, equal-length
#'   DNA sequences. Use uppercase A, T, C and G. For the supplied m6A model,
#'   provide five-nucleotide sequences (DNA 5-mers).
#'
#' @return A data frame with one row per input sequence and one factor column per
#'   nucleotide position, named `nt_pos1`, `nt_pos2`, and so on. Each column has
#'   levels `A`, `T`, `C`, `G`, in that order. Characters outside these levels
#'   become `NA`; the function does not validate sequence length or content.
#'
#' @examples
#' m6APrediction:::dna_encoding(c("AGACT", "TGACA"))
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A modification for multiple candidate sites
#'
#' Use a fitted random forest classifier to predict m6A modification probabilities
#' and statuses for the rows of a feature data frame. RNA categories are converted
#' to fixed factor levels and DNA 5-mers are encoded by nucleotide position before
#' prediction.
#'
#' @param ml_fit A fitted `randomForest` classification model trained on the
#'   corresponding features, including nucleotide-position factors. Its class
#'   probabilities must include a column named `Positive`.
#' @param feature_df A data frame with one row per candidate site and these
#'   required columns:
#'   \itemize{
#'     \item `gc_content`: numeric GC fraction, between 0 and 1.
#'     \item `RNA_type`: `mRNA`, `lincRNA`, `lncRNA`, or `pseudogene`.
#'     \item `RNA_region`: `CDS`, `intron`, `3'UTR`, or `5'UTR`.
#'     \item `exon_length`: numeric exon-length feature, already log2-transformed
#'       as in the course training data.
#'     \item `distance_to_junction`: numeric distance-to-junction feature, already
#'       log2-transformed as in the course training data.
#'     \item `evolutionary_conservation`: numeric conservation score on the same
#'       scale as the training data.
#'     \item `DNA_5mer`: a character column of uppercase, five-nucleotide DNA
#'       sequences containing only A, T, C and G.
#'   }
#'   Inputs should be non-missing and compatible with the fitted model.
#' @param positive_threshold A single numeric cutoff between 0 and 1. Defaults
#'   to 0.5. A site is `Positive` only when its predicted probability is strictly
#'   greater than this cutoff; equality is classified as `Negative`.
#'
#' @return A data frame with the original rows and columns plus
#'   `predicted_m6A_prob` (numeric probability for the `Positive` class) and
#'   `predicted_m6A_status` (character, `Positive` or `Negative`). The returned
#'   `RNA_type` and `RNA_region` columns are factors with the levels listed above.
#'   Existing prediction columns with these names are overwritten. Encoded
#'   nucleotide-position columns are used internally and are not added to the
#'   returned data frame.
#'
#' @details Only the presence of the required feature columns is explicitly
#'   checked. Unsupported RNA categories become `NA`; the caller must supply
#'   valid feature values. This function does not log-transform input features.
#'
#' @examples
#' library(m6APrediction)
#'
#' # Load the fitted model and example features bundled with the package.
#' rf_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' ))
#' feature_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv", package = "m6APrediction"
#' ))
#'
#' # Predict all candidate sites using a probability cutoff of 0.5.
#' predictions <- prediction_multiple(
#'   ml_fit = rf_fit,
#'   feature_df = feature_df,
#'   positive_threshold = 0.5
#' )
#' head(predictions[, c(
#'   "DNA_5mer", "predicted_m6A_prob", "predicted_m6A_status"
#' )])
#'
#' @seealso [prediction_single()]
#' @import randomForest
#' @importFrom stats predict
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict m6A modification for one candidate site
#'
#' Construct a one-row feature data frame and call [prediction_multiple()] to
#' obtain an m6A probability and status for a single site.
#'
#' @param ml_fit A fitted `randomForest` classification model compatible with
#'   these features and with a `Positive` class probability column.
#' @param gc_content A single numeric GC fraction between 0 and 1.
#' @param RNA_type A single RNA category: `mRNA`, `lincRNA`, `lncRNA`, or
#'   `pseudogene`.
#' @param RNA_region A single RNA region: `CDS`, `intron`, `3'UTR`, or `5'UTR`.
#' @param exon_length A single numeric exon-length feature, already
#'   log2-transformed as in the course training data.
#' @param distance_to_junction A single numeric distance-to-junction feature,
#'   already log2-transformed as in the course training data.
#' @param evolutionary_conservation A single numeric conservation score on the
#'   same scale as the training data.
#' @param DNA_5mer A single uppercase, five-nucleotide DNA sequence containing
#'   only A, T, C and G.
#' @param positive_threshold A single numeric cutoff between 0 and 1. Defaults
#'   to 0.5. Probabilities strictly greater than the cutoff are `Positive`;
#'   probabilities equal to or below it are `Negative`.
#'
#' @return A named character vector of length two, with elements
#'   `predicted_m6A_prob` and `predicted_m6A_status`. Combining the numeric
#'   probability and the status text coerces both elements to character. Use
#'   `as.numeric(result["predicted_m6A_prob"])` to recover a numeric probability.
#'
#' @details Supply one non-missing value per feature. Inputs must match the
#'   preprocessing used for model training; no log transformation is performed
#'   by this function.
#'
#' @examples
#' library(m6APrediction)
#'
#' # Load the fitted model and select one site from the example data.
#' rf_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"
#' ))
#' feature_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv", package = "m6APrediction"
#' ))
#' site <- feature_df[1, ]
#'
#' result <- prediction_single(
#'   ml_fit = rf_fit,
#'   gc_content = site$gc_content,
#'   RNA_type = site$RNA_type,
#'   RNA_region = site$RNA_region,
#'   exon_length = site$exon_length,
#'   distance_to_junction = site$distance_to_junction,
#'   evolutionary_conservation = site$evolutionary_conservation,
#'   DNA_5mer = site$DNA_5mer,
#'   positive_threshold = 0.5
#' )
#' result
#' # The single-site function returns text; convert the probability if needed.
#' as.numeric(result["predicted_m6A_prob"])
#'
#' @seealso [prediction_multiple()]
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
