# m6APrediction

## Introduction

This R package predicts m6A modification sites using a pre-trained Random Forest model. It was developed for the BIO215 practical course.

## Background

m6A (N6-methyladenosine) is an important RNA modification that regulates gene expression.

## Methods

- Feature preprocessing
- DNA sequence encoding
- Random Forest prediction
- Single-site and batch prediction

## Tools

- R
- RStudio
- randomForest

## Results

The package returns predicted m6A probabilities and Positive/Negative classifications. Model performance and evaluation results will be added here.

## Project Structure

```text
m6APrediction/
├── R/
│   └── main.R
├── inst/
│   └── extdata/
│       ├── rf_fit.rds
│       └── m6A_input_example.csv
├── man/
├── DESCRIPTION
├── NAMESPACE
└── README.md
```
