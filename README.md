# Extrapolating-Survival-Curves
Investigation of techniques for survival extrapolation using suitable richly parametric &amp; semi-parametric methods for time to event data. 

# Dataset
The dataset used in this analysis is sourced from [add dataset source here]. It contains clinical information about breast cancer patients, including demographic details, tumor characteristics, treatment history, and survival outcomes. The dataset comprises both numerical and categorical variables, providing a comprehensive view of patients' health status and disease progression.

# Code Structure
The code is written in R programming language and organized into several sections:

# Data Preprocessing: 
Reads the raw dataset, handles missing values, converts data types, and preprocesses variables for analysis.

# Exploratory Data Analysis (EDA): 
Visualizes the distribution of variables, explores relationships between predictors, and conducts statistical tests to identify associations using statistical tests such as Chi-Squared, Fisher's Exact Test, and Spearman Rank Correlation.

# Variable Selection: 
Several variable selection techniques are employed to identify the most relevant predictors for the survival analysis. This includes step-wise methods such as backward selection, forward selection, and both-direction selection.

# Principal Component Analysis (PCA): 
PCA is performed to reduce the dimensionality of the data and identify patterns among the predictors.

# Visualization: 
Generates visualizations including Kaplan-Meier plots, bar charts, histograms, scatter plots, box plots, and heatmaps to gain insights into the data and present findings.

# Modeling Techniques:
Models include Exponential, Weibull, Log-Logistic, Log-Normal, Gompertz, Gamma, and Generalized Gamma distributions.

# Model Evaluation and Selection:

- Computes and compares information criteria such as AIC and BIC for model selection.
Utilizes traceplots and autocorrelation plots to assess model convergence and diagnostics.
Advanced Modeling Techniques:

- Implements advanced models like Royston-Parmar splines and Poly-Weibull distributions.
Explores various model configurations and assesses their performance using AIC.
Extrapolation of Survival Curves:

- Extrapolates survival curves beyond the observed data range using different packages (e.g., rms and survival).

# README File: 
Provides an overview of the project, dataset, code structure, and instructions for running the analysis.
