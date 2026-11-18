# m6APrediction

A simple R package for predicting m6A modification using a trained random forest model. Developed for the BIO215 practical course.

## Installation

```r
install.packages("devtools")
devtools::install_github("tomaccc111/m6APrediction")
```

## Example

```r
library(m6APrediction)

model <- readRDS(system.file(
  "extdata", "rf_fit.rds", package = "m6APrediction"
))
features <- read.csv(system.file(
  "extdata", "m6A_input_example.csv", package = "m6APrediction"
))

results <- prediction_multiple(model, features)
head(results[, c("DNA_5mer", "predicted_m6A_prob", "predicted_m6A_status")])
```

## Functions

- `prediction_multiple()`: predict probabilities and statuses for multiple sites.
- `prediction_single()`: predict the probability and status for one site.

Input data requires `gc_content`, `RNA_type`, `RNA_region`, `exon_length`, `distance_to_junction`, `evolutionary_conservation`, and `DNA_5mer`. Length and distance features should already be log2-transformed, as in the example data.

The default probability threshold is 0.5: values above it are labelled `Positive`; all others are `Negative`. See `?prediction_multiple` and `?prediction_single` for details.
