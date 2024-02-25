# Extrapolating-Survival-Curves
Investigation of techniques for survival extrapolation using suitable richly parametric &amp; semi-parametric methods for time to event data. 

# Dataset
The dataset used in this analysis is sourced from Kaggle or other Machine Learning repositories (Breast Cancer METABRIC). It contains clinical information about breast cancer patients, including demographic details, tumor characteristics, treatment history, and survival outcomes. The dataset comprises both numerical and categorical variables, providing a comprehensive view of patients' health status and disease progression.

# Code Structure
The code is written in R programming language and organized into several sections:

# Data Preprocessing: 
Reads the raw dataset, handles missing values, converts data types, and preprocesses variables for analysis.

# Exploratory Data Analysis (EDA): 
Visualizes the distribution of variables, explores relationships between predictors, and conducts statistical tests to identify associations using statistical tests such as:
- Pearson Correlation,
- Chi-Squared,
- Fisher's Exact Test,
- Spearman Rank Correlation,
- Kendall's Tau Correlation,
- Mutual Information.

# Variable Selection: 
Several variable selection techniques are employed to identify the most relevant predictors for the survival analysis. This includes step-wise methods such as:
- backward selection,
- forward selection,
- both-direction selection,
- ridge regression,
- lasso regression,
- elastic net.

# Principal Component Analysis (PCA): 
PCA is performed to reduce the dimensionality of the data and identify patterns among the predictors.

# Visualization: 
Generates visualizations including Kaplan-Meier plots, bar charts, histograms, scatter plots, box plots, and heatmaps to gain insights into the data and present findings.

# Modeling Techniques:
Models fitted include both parametric & semi-parametric models, such as:
- Exponential,
- Weibull,
- Log-Logistic,
- Log-Normal,
- Gompertz,
- Gamma,
- and Generalized Gamma,
- Royston-Parmar splines,
- Poly-Weibull,
- Cox PH model.
Non-parametric estimators (Kaplan-Meier) are also implemented.

# Model Evaluation and Selection:

- Computes and compares information criteria such as AIC and BIC for model selection.

- Extrapolates survival curves beyond the observed data range using different packages (e.g., rms and survival).

- Calculates evaluation metrics such as the Mean Brier Score, L1 & L2 norms for each fitted model. These metrics help assess the goodness-of-fit and predictive performance of the models.

- Evaluate the calibration and discrimination of the fitted models.

- Calculate the mean survival time for each fitted model for different scenarios, both via closed-form solutions and numerical integration.


# README File: 
Provides an overview of the project, dataset, code structure, and instructions for running the analysis.
