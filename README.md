# Extrapolating-Survival-Curves
Investigation of techniques for survival extrapolation using suitable richly parametric &amp; semi-parametric methods for time to event data. 

MSc Thesis Project: Breast Cancer Analysis
Introduction
This repository contains the R script used for the analysis of breast cancer data as part of my MSc thesis project. The project focuses on exploring various clinical features and their association with overall survival in breast cancer patients.

Dataset
The dataset used in this analysis is the "brca_metabric_clinical_data_1980.xlsx" file, which contains clinical data of 1980 breast cancer patients. The dataset includes 39 variables, including patient demographics, tumor characteristics, treatment details, and overall survival status.

Analysis
Data Preprocessing
Loaded the dataset using the readxl library.
Checked for missing values and converted character missing values to NA.
Omitted rows with missing values.
Variable Selection
Created variables of interest, including overall survival status, age at diagnosis, type of breast surgery, cancer type, cancer cellularity, chemotherapy status, hormonal therapy status, etc.
Exploratory Data Analysis (EDA)
Visualized missing values using the naniar and ggplot2 libraries.
Visualized categorical variables using bar plots and contingency tables.
Visualized numerical variables using histograms and scatter plots.
Statistical Analysis
Conducted Chi-squared tests and Fisher's exact tests to examine associations between categorical variables.
Calculated Pearson correlation coefficients to examine correlations between numerical variables.
