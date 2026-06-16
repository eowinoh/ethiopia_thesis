# Maternal Dietary Diversity and Child Nutritional Outcomes Analysis

This repository contains the code, statistical analyses, and supporting materials used to investigate factors associated with maternal dietary diversity and child nutritional outcomes. The analyses examine both maternal and child nutrition using a range of regression approaches, including Poisson regression with robust standard errors, Conway-Maxwell-Poisson (COM-Poisson) regression, and quantile regression.

## Objectives

The study aims to:

* Identify factors associated with maternal dietary diversity.
* Examine determinants of child nutritional outcomes, including:

  * Length-for-age z-score (LAZ)
  * Weight-for-age z-score (WAZ)
  * Weight-for-length z-score (WLZ)
* Assess whether associations vary across different points of the outcome distributions using quantile regression.


## Methods

The repository includes code for:

* Data cleaning and preprocessing
* Missing data handling using multiple imputation
* Descriptive statistical analyses
* Poisson regression with robust standard errors
* Conway-Maxwell-Poisson (COM-Poisson) regression
* Quantile regression

## Order of Execution
```text
Data_Preprocessing.R
        ↓
Final_Analysis_Processing.R
        ↓
  Objective_1.R
        ↓
  Objective_2.R
        ↓
  Objective_3.R       
```
## Software

Analyses were conducted in R using packages for count-data modelling, multiple imputation, and quantile regression.

## Reproducibility

All analyses are designed to be reproducible from the provided scripts. Users should ensure that the required packages are installed and that data access permissions are respected where datasets are not publicly available.

## License

This repository is provided for academic and research purposes.

```
```
