

library(readxl)
data <- read_excel("C:/Users/mihal/OneDrive/brca_metabric_clinical_data_1980.xlsx",col_names = T) # 1980 obs of 39 variables
head(data)
dim(data) # 1980 * 39
str(data) # 21 Clinical Features
length(which(is.na(data))) # No values coded as NA, however there are missing values!

# Convert the character missing values
data[ data == "NA" ] <- NA

# Data left after omitting the missing values
df <- na.omit(data) 
str(df) # 1092 obs

# Creating our variables
event <- df$`Overall Survival Status` # Target variable whether the patient is alive of dead.
table(event)

event[event == '0:LIVING'] <- 0
event[event == '1:DECEASED'] <- 1
event <- as.numeric(event)
table(event)

age <- df$`Age at Diagnosis` # Age of the patient at diagnosis time

type_breast_surg <- factor(df$`Type of Breast Surgery`) # Breast cancer surgery type: 1- MASTECTOMY, which refers to a surgery to 
# remove all breast tissue from a breast as a way to treat or prevent breast cancer. 2- BREAST CONSERVING, which refers to a urgery 
# where only the part of the breast that has cancer is removed

cancer_type <- factor(df$`Cancer Type`) # Breast cancer types: 1- Breast Cancer or 2- Breast Sarcoma
table(cancer_type) # 1 level!

cancer_type_det <- factor(df$`Cancer Type Detailed`)

cell <- factor(df$Cellularity) # Cancer cellularity post chemotherapy, which refers to the amount of tumor cells in the specimen and 
# their arrangement into clusters

chemo <- factor(df$Chemotherapy) # Whether or not the patient had chemotherapy as a treatment (yes/no)

Pam50 <- factor(df$`Pam50 + Claudin-low subtype`) # Pam 50: is a tumor profiling test that helps show whether some estrogen 
# receptor-positive (ER-positive), HER2-negative breast cancers are likely to metastasize (when breast cancer spreads to other organs). 
# The claudin-low breast cancer subtype is defined by gene expression characteristics, most prominently: Low expression of cell–cell
# adhesion genes, high expression of epithelial–mesenchymal transition (EMT) genes, and stem cell-like/less differentiated gene 
# expression patterns

cohort <- factor(df$Cohort) # Cohort is a group of subjects who share a defining characteristic (It takes a value from 1 to 5)

er_ihc <- factor(df$`ER status measured by IHC`) # To assess if estrogen receptors are expressed on cancer cells by using
# immune-histochemistry (a dye used in pathology that targets specific antigen, if it is there, it will give a color, it is not there, 
# the tissue on the slide will be colored) (positive/negative)

er <- factor(df$`ER Status`) # Cancer cells are positive or negative for estrogen receptors

nhg <- factor(df$`Neoplasm Histologic Grade`) # Determined by pathology by looking the nature of the cells, do they look aggressive or
# not (It takes a value from 1 to 3)

her2_snp6 <- factor(df$`HER2 status measured by SNP6`) # To assess if the cancer positive for HER2 or not by using advance molecular 
# techniques (Type of next generation sequencing)

her2 <- factor(df$`HER2 Status`) # Whether the cancer is positive or negative for HER2

tohs <- factor(df$`Tumor Other Histologic Subtype`) # Type of the cancer based on microscopic examination of the cancer tissue 
# (It takes a value of 'Ductal/NST', 'Mixed', 'Lobular', 'Tubular/ cribriform', 'Mucinous', 'Medullary', 'Other', 'Metaplastic' )

hormone <- factor(df$`Hormone Therapy`) # Whether or not the patient had hormonal as a treatment (yes/no)

ims <- factor(df$`Inferred Menopausal State`) # Whether the patient is is post menopausal or not (post/pre)

int_clust <- factor(df$`Integrative Cluster`) # Molecular subtype of the cancer based on some gene expression 
# (It takes a value from '4ER+', '3', '9', '7', '4ER-', '5', '8', '10', '1', '2', '6')

ptl <- factor(df$`Primary Tumor Laterality`) # Whether it is involving the right breast or the left breast

lnep <- as.numeric(df$`Lymph nodes examined positive`) # To take samples of the lymph node during the surgery and see if there were 
# involved by the cancer

mc <- as.numeric(df$`Mutation Count`) # Number of gene that has relevant mutations

npi <- df$`Nottingham prognostic index` # It is used to determine prognosis following surgery for breast cancer. Its value is 
# calculated using three pathological criteria: the size of the tumour; the number of involved lymph nodes; and the grade of the tumour.

oc <- factor(df$`Oncotree Code`) # The OncoTree is an open-source ontology that was developed at Memorial Sloan Kettering Cancer 
# Center (MSK) for standardizing cancer type diagnosis from a clinical perspective by assigning each diagnosis a unique OncoTree code.

Time <- df$`Overall Survival (Months)` # Duration from the time of the intervention to death

prs <- factor(df$`PR Status`) # Cancer cells are positive or negative for progesterone receptors

radio <- factor(df$`Radio Therapy`) # Whether or not the patient had radio as a treatment (yes/no)

sex <- factor(df$Sex) # Careful there's 1 level

gcs <- factor(df$`3-Gene classifier subtype`) # Three Gene classifier subtype It takes a value from 'ER-/HER2-', 
# 'ER+/HER2- High Prolif', nan, 'ER+/HER2- Low Prolif','HER2+'

tmb <- df$`TMB (nonsynonymous)` # TMB, defined as the number of somatic mutations per megabase of interrogated genomic sequence,
# varies across malignancies

tumor_size <- as.numeric(df$`Tumor Size`) # Tumor size measured by imaging techniques

tumor_stage <- factor(df$`Tumor Stage`) # Stage of the cancer based on the involvement of surrounding structures, lymph nodes and 
# distant spread

#vital_status <- df$`Patient's Vital Status` # Died of disease etc




# Note: Cancer grade vs tumor stage

# What is cancer grade?
# A cancer’s grade describes how abnormal the cancer cells and tissue look under a microscope when compared to healthy cells. 
# Cancer cells that look and organize most like healthy cells and tissue are low grade tumors. Doctors describe these cancers as being
# well differentiated. Lower grade cancers are typically less aggressive and have a better prognosis.
# The more abnormal the cells look and organize themselves, the higher the cancer’s grade. Cancer cells with a high grades tend to be 
# more aggressive. They are called poorly differentiated or undifferentiated.
# Some cancers have their own system for grading tumors. Many others use a standard 1-4 grading scale.
# Grade 1: Tumor cells and tissue looks most like healthy cells and tissue. These are called well-differentiated tumors and are 
# considered low grade.
# Grade 2: The cells and tissue are somewhat abnormal and are called moderately differentiated. These are intermediate grade tumors.
# Grade 3: Cancer cells and tissue look very abnormal. These cancers are considered poorly differentiated, since they no longer have an
# architectural structure or pattern. Grade 3 tumors are considered high grade.
# Grade 4: These undifferentiated cancers have the most abnormal looking cells. These are the highest grade and typically grow and 
# spread faster than lower grade tumors.

# What is a cancer stage?
# While a grade describes the appearance of cancer cells and tissue, a cancer’s stage explains how large the 
# primary tumor is and how far the cancer has spread in the patient’s body.
# There are several different staging systems. Many of these have been created for specific kinds of cancers. Others can be used to 
# describe several types of cancer.
# Stage 0 to stage IV
# One common system that many people are aware of puts cancer on a scale of 0 to IV.
# Stage 0 is for abnormal cells that haven’t spread and are not considered cancer, though they could become cancerous in the future. 
# This stage is also called “in-situ.”
# Stage I through Stage III are for cancers that haven’t spread beyond the primary tumor site or have only spread to nearby tissue. The
# higher the stage number, the larger the tumor and the more it has spread.
# Stage IV cancer has spread to distant areas of the body.






##############################################       Visualize the NAs


library(naniar)
# Percentage of the total missing values in the data set
pct_miss(data) # 1.71 %


library(naniar)
library(ggplot2)

gg_miss_var(data) +
  theme(
    panel.background = element_rect(fill = "grey"),
    text = element_text(color = "black")
  )


# Data visualizations using package naniar
gg_miss_var(data) +
  theme(text = element_text(color = "black"))


miss_var_summary(data, order = T)


# Provide a summary for each variable of the number, percent missings, and cumulative sum of missings of the order of the variables. 
# By default, it orders by the most missings in each variable.

library(ggplot2)
# Where are missings located?
vis_miss(data) + theme(axis.text.x = element_text(angle=80))
 


# Change the color of the fill
#  vis_miss(data) + 
#  scale_fill_manual(values = c("grey", "black")) + 
#  scale_color_manual(values = c("black", "black")) +
#  theme(axis.text.x = element_text(angle=80))






##############################################       Merging some levels


table(tohs)

# Merge tohs levels
tohs[tohs == 'Medullary'] <- 'Other'
tohs[tohs == 'Mucinous'] <- 'Other'
tohs[tohs == 'Tubular/ cribriform'] <- 'Other'
table(tohs)

tohs <- factor(tohs,levels = c('Ductal/NST','Lobular','Mixed','Other'))
table(tohs)

# Merge oc levels
levels(oc) <- c(levels(oc),'Other')
oc[oc == 'BREAST'] <- 'Other'
oc[oc == 'IMMC'] <- 'Other'
oc[oc == 'ILC'] <- 'Other'
table(oc)

oc <- factor(oc, levels=c('IDC','MDLC','Other'))
table(oc)




# Merge cancer_type_det levels
levels(cancer_type_det) <- c(levels(cancer_type_det),'Other')
cancer_type_det[cancer_type_det == 'Breast'] <- 'Other'
cancer_type_det[cancer_type_det == 'Breast Invasive Mixed Mucinous Carcinoma'] <- 'Other'
cancer_type_det[cancer_type_det == 'Invasive Mixed Mucinous Carcinoma'] <- 'Other'

table(cancer_type_det)

cancer_type_det <- factor(cancer_type_det, levels=c('Breast Invasive Ductal Carcinoma',
                                                    'Breast Invasive Lobular Carcinoma','Breast Mixed Ductal and Lobular Carcinoma',
                                                    'Other'))

table(cancer_type_det)





# Create a data frame with all the predictors 
predictors <- data.frame( age = age, type_breast_surg = type_breast_surg,
cancer_type_det = cancer_type_det,  
cell = cell, chemo = chemo, Pam50 = Pam50, cohort = cohort, er = er,  nhg = nhg, her2 = her2, tohs = tohs, 
hormone = hormone, ims = ims, ptl = ptl, lnep = lnep, mc = mc, npi = npi, oc = oc,
prs = prs, radio = radio, gcs = gcs, tumor_size = tumor_size, tumor_stage = tumor_stage, tmb = tmb,
Time = Time, event = event )

str(predictors) 
# 26 predictors (including Time and Event) without: tmb (correlated with mc), int_clust, her2_snp6, er_ihc, study ID, patient ID,
# relapse free status x2, Number of Samples Per Patient, Sample Type, Sex (1 level), Patient's Vital Status, cancer_type (1 level)







##############################################       Cut off data along with the desirable KM plot


# number of patients
n_pts <- dim(predictors)[1]
n_pts

# new dataset: observation time <= 160 months for all subjects 
# (therefore for all pts with Time > 160 we must have event = 0, everything else stays the same)
cut_off_dat <- predictors
for (i in 1:n_pts){
  if (predictors$Time[i] > 160){ 
    cut_off_dat$Time[i] <- 160
    cut_off_dat$event[i] <- 0
  }
}

# check
sum(cut_off_dat$Time > 160) # 0 cases -> ok!
sum(predictors$Time > 160) # initially, 367 cases had Time > 160 months
sum(cut_off_dat$Time == 160) # also 368 cases as expected
# we now check if the event indicator changed 
x <- predictors[predictors$Time > 160,]
sum(x$event) # 115 events after 160 months
length(which(cut_off_dat$event != predictors$event)) # indeed we have 115 cases with different indicator 

library(survival)
library(ggfortify)

km <- survfit(Surv(Time, event) ~ 1, type = "kaplan-meier", data = cut_off_dat)

plot(km , xlab="Months", ylab='S(t)', conf.int = F, main = 'Kaplan Meyer Plot',xlim = c(0,337), lwd = 2, las = 1)
# autoplot(km ,conf.int = T, xlab = 'Months', main = 'Kaplan Meier Plot', xlim=c(0,337) )
# print(km,print.rmean = T)

abline(h = 0.5, lty = 2)

abline(v = 160, lty = 2, col = 'red')

str(cut_off_dat) # 26 predictors


# Export cut_off_dat to CSV file with the name 'predictors'
# write.csv(cut_off_dat,file='predictors.csv')






##############################################       Interactive KM plot via plotly


event.full.data <- data$`Overall Survival Status` 
table(event.full.data)

event.full.data[event.full.data == '0:LIVING'] <- 0
event.full.data[event.full.data == '1:DECEASED'] <- 1
event.full.data <- as.numeric(event.full.data)
table(event.full.data)

time <- data$`Overall Survival (Months)`

library(survminer)
library(plotly)



km.full <- survfit(Surv(time, event.full.data) ~ 1, type = "kaplan-meier")

# Interactive KM curve on the full data
ppp <- ggsurvplot(km.full, data = data)
ggplotly(ppp[[1]])



##############################################       Visualization of our variables



# Creating a data frame with the categorical variables
df_categorical_preds <- data.frame(type_breast_surg = cut_off_dat$type_breast_surg, cancer_type_det = cut_off_dat$cancer_type_det,
                                   cell = cut_off_dat$cell, chemo = cut_off_dat$chemo, Pam50  = cut_off_dat$Pam50, cohort = cut_off_dat$cohort,
                                   er = cut_off_dat$er, nhg = cut_off_dat$nhg, her2 = cut_off_dat$her2, tohs = cut_off_dat$tohs,
                                   hormone = cut_off_dat$hormone, ims = cut_off_dat$ims, ptl = cut_off_dat$ptl,
                                   oc = cut_off_dat$oc, prs = cut_off_dat$prs, radio = cut_off_dat$radio, 
                                   gcs = cut_off_dat$gcs, tumor_stage = cut_off_dat$tumor_stage)



# Creating a data frame with the numerical variables
df_numerical_preds <- data.frame( age = age, tumor_size = tumor_size, lnep = lnep, mc = mc, tmb = tmb, npi = npi )






# Function to plot bar charts of the categorical variables
stacked_bar_charts <- function(data) {
  # loop through each variable in the data frame
  for (col in colnames(data)) {
    # create a table of counts for each level of the variable
    counts <- table(data[, col])
    # create a stacked bar chart of the counts
    barplot(counts, main = col, xlab = col, col = rainbow(length(counts)))
  }
}

stacked_bar_charts(df_categorical_preds)






library(ggplot2)

# Function to plot histograms of numerical variables
coloured_histograms <- function(df) {
  # Identify numerical variables in the data frame
  numerical_vars <- sapply(df, is.numeric)
  
  # Create a list of the numerical variable names
  var_names <- names(df)[numerical_vars]
  
  # Set the number of columns for the plot grid
  num_cols <- ceiling(sqrt(length(var_names)))
  
  # Set up the plot grid
  par(mfrow = c(num_cols, num_cols))
  
  # Loop through each numerical variable
  for (i in 1:length(var_names)) {
    # Create a histogram with 50 bins and a coloured fill
    hist(df[, var_names[i]], breaks = 50, col = "cornflowerblue", main = var_names[i])
  }
}

coloured_histograms(df_numerical_preds)





# Function to create scatter plots for each pair of numerical variables
scatter_plot_pairs <- function(df) {
  num_cols <- ncol(df)
  
  # Loop over each pair of variables
  for (i in 1:(num_cols - 1)) {
    for (j in (i + 1):num_cols) {
      # Create the scatter plot for this pair of variables
      plot(df[,i], df[,j], xlab=colnames(df)[i], ylab=colnames(df)[j])
    }
  }
}

scatter_plot_pairs(df_numerical_preds)








# Function to plot box plots of each variable of the data set
plot_boxplots <- function(data) {
  # Loop through each variable in the data frame
  for (col in names(data)) {
    # Create a boxplot for the variable
    boxplot(data[[col]], main = col)
  }
}


plot_boxplots(df_numerical_preds)



# Function to plot box plots of each variable of the data set in the same grid 
coloured_boxplots <- function(df) {
  # Identify numerical variables in the data frame
  numerical_vars <- sapply(df, is.numeric)
  
  # Create a list of the numerical variable names
  var_names <- names(df)[numerical_vars]
  
  # Set the number of columns for the plot grid
  num_cols <- ceiling(sqrt(length(var_names)))
  
  # Set up the plot grid
  par(mfrow = c(num_cols, num_cols))
  
  # Loop through each numerical variable
  for (i in 1:length(var_names)) {
    # Create a box plot with a coloured fill
    boxplot(df[, var_names[i]], col = "cornflowerblue", main = var_names[i])
  }
}


coloured_boxplots(df_numerical_preds)







##############################################       Contingency tables, X^2 and Fischer tests


# Contingency tables and Chi-Squared tests
tab1 <- table(her2_snp6,her2)
tab1

library(stats)
fisher.test(tab1)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 



tab2 <- table(cancer_type_det,type_breast_surg)
tab2
prop.test(tab2) # p-value = 0.2011 => fail to reject H0 => the two categorical variables are independent 
# There's no association between the two variables



tab3 <- table(Pam50,her2)
tab3

prop.test(tab3)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 

chisq.test(tab3)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 



tab4 <- table(Pam50, her2_snp6)
tab4

chisq.test(tab4)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 


tab5 <- table(er_ihc,er)
tab5

prop.test(tab5) 
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 

fisher.test(tab5)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 



tab6 <- table(gcs,er)
tab6

prop.test(tab6)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 

chisq.test(tab6)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 



tab7 <- table(gcs,er_ihc)
tab7

prop.test(tab7)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 




tab8 <- table(Pam50, gcs)
tab8

chisq.test(tab8)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 



tab9 <- table(Pam50, nhg)
tab9 

chisq.test(tab9)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 



tab10 <- table(gcs, nhg)
tab10

chisq.test(tab10)
# p_value ~= 0 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables). 



tab11 <- table(Pam50,type_breast_surg)
tab11 

chisq.test(tab11)
# p_value = 0.03 => Reject H0 => The two categorical variables are NOT independent (there's association between the two variables).



tab12 <- table(type_breast_surg, tohs)
tab12


chisq.test(tab12)
# p_value = 0.11 => Fail to reject H0 => The two categorical variables are independent (there's no association between the two variables).




tab13 <- table(type_breast_surg, oc)
tab13


chisq.test(tab13)
# p_value = 0.28 => Fail to reject H0 => The two categorical variables are independent (there's no association between the two variables).



#install.packages("pander")
library(pander)

# Create a pander table from the matrix
pander_table <- pandoc.table(res, style = "rmarkdown")

# Print the pander table
print(pander_table)







library(gplots)

# Function to visualize the cells that have p.value <= 0.05 and add a title
create_heatmap <- function(data, title) {
  # Set the color scale for the heatmap
  my_col <- colorRampPalette(c("white", "blue"))(256)
  
  # Create a new plot window and set the margins
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # Create the heatmap using the heatmap.2 function
  heatmap.2(data,
            dendrogram = "none",
            scale = "none",
            trace = "none",
            col = my_col,
            key = TRUE,
            keysize = 1.0,
            symkey = FALSE,
            density.info = "none",
            Rowv = FALSE,
            Colv = FALSE,
            cexCol = 1,
            cexRow = 1,
            cellnote = data <= 0.05,  # Highlight cells where value is below 0.05
            notecol = "black",
            notecex = 1,
            key.title = NA)
  
  # Add a title to the plot
  title(main = title)
}

# Call the modified function with the desired title
create_heatmap(res, "Chi-Squared Test of Independence")








# Function to visualize the cells that have p.value <= 0.05 and add a title
create_heatmap <- function(data, title, fontsize = 12) {
  # Set the color scale for the heatmap
  my_col <- colorRampPalette(c("white", "blue"))(256)
  
  # Create a new plot window and set the margins
  par(mar = c(5, 5, 4, 2) + 0.1)
  
  # Create the heatmap using the heatmap.2 function
  heatmap.2(data,
            dendrogram = "none",
            scale = "none",
            trace = "none",
            col = my_col,
            key = TRUE,
            keysize = 1.0,
            symkey = FALSE,
            density.info = "none",
            Rowv = FALSE,
            Colv = FALSE,
            cexCol = fontsize / 10,  # Adjust column label font size
            cexRow = fontsize / 10,  # Adjust row label font size
            cellnote = data <= 0.05,  # Highlight cells where value is below 0.05
            notecol = "black",
            notecex = fontsize / 10,
            key.title = NA)
  
  # Add a title to the plot
  title(main = title, cex.main = fontsize / 10)
}

# Call the modified function with the desired title and font size
create_heatmap(res, "Chi-Squared Test of Independence", fontsize = 14)









##############################################       Pearson Correlation for the numerical variables



library(corrplot)
## corrplot 0.92 loaded
new.corrmat <- cor(as.matrix(df_numerical_preds))
corrplot(new.corrmat, method = 'number', tl.cex = 1.5)




##############################################       Mutual Information

# The Mutual Information between two random variables measures non-linear relations between them. Besides, it indicates how much 
# information can be obtained from a random variable by observing another random variable.
# It is closely linked to the concept of entropy. This is because it can also be known as the reduction of uncertainty of a random 
# variable if another is known. Therefore, a high mutual information value indicates a large reduction of uncertainty whereas a low 
# value indicates a small reduction. If the mutual information is zero, that means that the two random variables are independent.

# The main difference is that correlation is a measure of linear dependence, whereas mutual information measures general dependence 
# (including non-linear relations). Therefore, mutual information detects dependencies that do not only depend on the covariance.
# Thus, mutual information is zero when the two random variables are strictly independent.


library(infotheo)


compute_mi <- function(data) {
  # get the names of all variables in the data frame
  variable_names <- names(data)
  
  # initialize an empty matrix to store the correlations
  mi_matrix <- matrix(NA, nrow = length(variable_names), ncol = length(variable_names))
  
  # loop through each pair of variables
  for (i in 1:length(variable_names)) {
    for (j in 1:length(variable_names)) {
      # don't compute the correlation between a variable and itself
      if (i == j) {
        mi_matrix[i, j] <- 1
      } else {
        # Calculate mutual information between the two variables (the default is the "emp" : This estimator computes the entropy of the empirical probability distribution.)
        mi_matrix[i, j] <- mutinformation(data[[i]], data[[j]])
      }
    }
  }
  
  # name the rows and columns of the correlation matrix with the variable names
  rownames(mi_matrix) <- variable_names
  colnames(mi_matrix) <- variable_names
  
  # return the correlation matrix
  return(mi_matrix)
}


mi_mat <- round(compute_mi(df_categorical_preds[,-19]), 3)
mi_mat



# Heatmap of the correlations of our dataset
library(corrplot)

#png(filename = "mycorrplot.png", width = 1500, height = 1500)

corrplot(mi_mat, method = "number", addCoef.col = 1,number.cex = 1)







##############################################      Spearman Rank Correlation


# Function to calculate Spearman Rank correlation between pairs of numerical variables
calculate_spearman_cor <- function(data) {
  
  # Create empty data frame to store results
  result_df <- data.frame(variable1 = character(),
                          variable2 = character(),
                          correlation = numeric(),
                          stringsAsFactors = FALSE)
  
  # Get list of numerical variables in data frame
  numeric_vars <- sapply(data, is.numeric)
  
  # Loop over each pair of numerical variables and calculate correlation
  for (i in 1:(length(numeric_vars) - 1)) {
    for (j in (i + 1):length(numeric_vars)) {
      var1 <- names(numeric_vars)[i]
      var2 <- names(numeric_vars)[j]
      cor_value <- cor(data[[var1]], data[[var2]], method = "spearman")
      
      # Add result to data frame
      result_df <- rbind(result_df, data.frame(variable1 = var1,
                                               variable2 = var2,
                                               correlation = cor_value,
                                               stringsAsFactors = FALSE))
    }
  }
  
  # Sort results by correlation value in descending order
  result_df <- result_df[order(-result_df$correlation),]
  
  # Print results
  print(result_df)
}


df_spearman_1 <- calculate_spearman_cor(df_numerical_preds)
df_spearman_1



# Create the data frame
df_spearman_final <- data.frame(Variable1 = c("mc", "lnep", "tumor_size", "tumor_size", "age"),
                 Variable2 = c("tmb", "npi", "npi", "lnep", "tumor_size"),
                 Correlation = c(0.998501316, 0.772848211, 0.594195026, 0.359519997, 0.123968171))
df_spearman_final




library(kableExtra)

# Create a fancy matrix with the mean survival time estimates for the 50% scenario
caption.5 <- "<b>Top 5 Spearman rank correlations in descending order</b>"
# Create a styled table using kableExtra
table.5 <- kable(df_spearman_final, format = "html", table.attr = "class='table'", 
                 caption = caption.5) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Display the table
table.5






##############################################      Kendall Rank Correlation



# Function to create contingency tables for each pair of variables, test for Kendall's correlation and store the p.values
  kendall_rank_test <- function(data) {
  n_col <- ncol(data)
  output_data <- data.frame(matrix(ncol = n_col, nrow = n_col))
  
  for (i in 1:n_col) {
    for (j in 1:n_col) {
      if (i == j) {
        output_data[i, j] <- NA
      } else {
        test <- cor.test(data[,i], data[,j], method = 'kendall')
        output_data[i, j] <- test$p.value
      }
    }
  }
  
  colnames(output_data) <- colnames(data)
  rownames(output_data) <- colnames(data)
  return(output_data)
}

 

  

# Round and store the results in a data frame named res
results_kendall <- round(kendall_rank_test(df_numerical_preds),3)
results_kendall




library(kableExtra)

# Create a fancy matrix with the mean survival time estimates for the 50% scenario
caption.6 <- "<b>Kendall's tau correlation</b>"
# Create a styled table using kableExtra
table.6 <- kable(results_kendall, format = "html", table.attr = "class='table'", 
                 caption = caption.6) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)
 

# Display the table
table.6




##############################################      Step-wise methods for variable selection


library(survival)
mod.cox.1 <- coxph(formula = Surv(Time,event) ~ age + type_breast_surg + cancer_type_det + tumor_size
+ cell + chemo + Pam50 + cohort + er + nhg + her2 + tohs + hormone + ims + ptl + lnep + mc + npi + oc + prs + radio
+ gcs + tumor_stage ) # 25 predictors (including time and event)


summary(mod.cox.1)
# tohs, oc -> NA's


# Fitting a Cox model after removing tohs and oc
mod.cox.2 <- coxph(formula = Surv(Time,event) ~ age + type_breast_surg + cancer_type_det + tumor_size
+ cell + chemo + Pam50 + cohort + er + nhg + her2 + hormone + ims + ptl + lnep + mc + npi + prs + radio 
+ gcs + tumor_stage)

summary(mod.cox.2) # 23 predictors (including Time and Event) - without -> tohs and oc



# Stepwise model selection
bth <- step(mod.cox.2, direction = "both",  trace = 0) # 11 predictors!
bth$coefficients


# Backward model selection
bckwrd <- step(mod.cox.2, direction = "backward") # 11 predictors!


# The 11 most important predictors from the 2 above step-wise procedures
formula = Surv(Time, event) ~ age + cancer_type_det + tumor_size + lnep + nhg + her2 + ims + npi + radio + gcs + hormone



# Create a new df with the 23 predictors
dat <- data.frame( age = age, type_breast_surg = type_breast_surg,
                   cancer_type_det = cancer_type_det,  
                   cell = cell, chemo = chemo, Pam50 = Pam50, cohort = cohort, er = er,  nhg = nhg, her2 = her2,
                   hormone = hormone, ims = ims, ptl = ptl, lnep = lnep, mc = mc, npi = npi,
                   prs = prs, radio = radio, gcs = gcs, tumor_size = tumor_size, tumor_stage = tumor_stage, tmb = tmb)


str(dat) # 22 predictors



# Forward model selection
cox.0 <- coxph(formula = Surv(Time,event) ~ 1, data = dat)
step(cox.0, direction = "forward", scope = list(lower = formula(cox.0),upper = formula(mod.cox.2))) # 10 predictors! -> cohort is the difference

# The 10 most important predictors from the forward model selection
formula = Surv(Time, event) ~ age + npi + cancer_type_det + tumor_size + lnep + cohort + her2 + hormone + ims + chemo




age_star = c("*","*","*")
type_breast_surg_star = c("","","")
cancer_type_det_star= c("*","*","*")
cell_star = c("","","")
chemo_star = c('',"","*")
Pam50_star = c('',"","")
cohort_star = c('',"","")
er_star = c('',"","")
nhg_star = c("*","*","")
her2_star = c("*","*","*")
hormone_star = c("*","*","*")
ims_star = c("*","*","*")
ptl_star = c('',"","")
lnep_star = c("*","*","*")
mc_star = c('',"","")
npi_star = c("*","*","*")
prs_star = c('',"","")
radio_star = c("*","*","*")
gcs_star = c("*","*","")
tumor_size_star = c("*","*","*")
tumor_stage_star = c("","","")
tmb_star = c('',"","")

var_sel_mat <- matrix( c(age_star, type_breast_surg_star,
                      cancer_type_det_star,  
                      cell_star, chemo_star, Pam50_star, cohort_star, er_star, nhg_star, her2_star,
                      hormone_star, ims_star, ptl_star, lnep_star, mc_star, npi_star,
                      prs_star, radio_star, gcs_star, tumor_size_star, tumor_stage_star, tmb_star), 
                      nrow = 3, ncol = 22, byrow = F
)

var_sel_mat

is.matrix(var_sel_mat)


colnames(var_sel_mat) <- c("Age", "Type_breast_surg","Cancer_type_det","Cell","Chemo","Pam50",
                       "Cohort","Er","Nhg","Her2","Hormone","Ims","Ptl","Lnep","Mc","Npi","Prs","Radio","Gcs","Tumor_size","Tumor_stage","Tmb"                   
                                            
                                          )
rownames(var_sel_mat) <- c("Backward","Both","Forward")
var_sel_mat




library(kableExtra)

# Convert matrix to a data frame
var_sel_df <- as.data.frame(var_sel_mat)




# Create a fancy matrix
caption.2.1 <- "<b>Stepwise Variable Selection</b>"
# Create a styled table using kableExtra
table.2.1 <- kable(var_sel_df, format = "html", table.attr = "class='table'", 
                 caption = caption.2.1) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Display the table
table.2.1





##############################################      Vif check for the variables selected by the step-wise methods


mod1 <- lm(Time ~ age + cancer_type_det + tumor_size + lnep + nhg + her2 + hormone + ims + npi + radio + gcs) # 11 predictors
summary(mod1)


library(car)
# Create vector of VIF values
vif_values <- vif(mod1)
vif_values


# We can also check for the presence of multicollinearity by looking at the tolerance values. The tolerance is simply the inverse of 
# the VIF, and a tolerance value less than 0.1 indicates a high degree of multicollinearity.
tolerance <- 1/vif_values
tolerance


# Corrected GVIF
GVIF_sqrt <- sqrt(vif_values[,3])
GVIF_sqrt # No value > 2



library(rms) # Regression Modeling Strategies

d <- datadist(predictors)
options(datadist = "d")


# Fitting a Cox PH model with the covariates selected by the automatic procedures (not all!)
cox.mod.rms <- cph(Surv(Time, event) ~  age + cancer_type_det + tumor_size + lnep + nhg + her2 + hormone + ims + npi + radio + gcs, surv=T)



cox_vif <- vif(cox.mod.rms)
cox_vif




# Sort the numerical values from largest to smallest
sorted_values <- sort(GVIF_sqrt, decreasing = TRUE)


# Sort the numerical values from largest to smallest
sorted_values <- sort(vif_values, decreasing = TRUE)


# Create a bar plot with sorted values
barplot(sorted_values, col = "skyblue", main = "Sorted VIF values",
        xlab = "Predictors", ylab = "VIF",  space = 0.5, las = 1,  cex.axis = 1.3)


##############################################      Dummy coding


cut_off_dat$Pam50_1 <- ifelse(cut_off_dat$Pam50 == 'Basal', 1, 0)
cut_off_dat$Pam50_2 <- ifelse(cut_off_dat$Pam50 == 'claudin-low', 1, 0)
cut_off_dat$Pam50_3 <- ifelse(cut_off_dat$Pam50 == 'Her2', 1, 0)
cut_off_dat$Pam50_4 <- ifelse(cut_off_dat$Pam50 == 'LumA', 1, 0)
cut_off_dat$Pam50_5 <- ifelse(cut_off_dat$Pam50 == 'LumB', 1, 0)
cut_off_dat$Pam50_6 <- ifelse(cut_off_dat$Pam50 == 'Normal', 1, 0)
#cut_off_dat$Pam50_7 <- ifelse(cut_off_dat$Pam50 == 'NC', 1, 0)

Pam50_1 <- cut_off_dat$Pam50_1
Pam50_2 <- cut_off_dat$Pam50_2
Pam50_3 <- cut_off_dat$Pam50_3
Pam50_4 <- cut_off_dat$Pam50_4
Pam50_5 <- cut_off_dat$Pam50_5
Pam50_6 <- cut_off_dat$Pam50_6
#Pam50_7 <- cut_off_dat$Pam50_7


#cut_off_dat$cancer_type_det_1 <- ifelse(cut_off_dat$cancer_type_det == 'Breast Invasive Ductal Carcinoma', 1, 0)
cut_off_dat$cancer_type_det_2 <- ifelse(cut_off_dat$cancer_type_det == 'Breast Invasive Lobular Carcinoma', 1, 0)
cut_off_dat$cancer_type_det_3 <- ifelse(cut_off_dat$cancer_type_det == 'Breast Mixed Ductal and Lobular Carcinoma', 1, 0)
cut_off_dat$cancer_type_det_4 <- ifelse(cut_off_dat$cancer_type_det == 'Other', 1, 0)

#cancer_type_det_1 <- cut_off_dat$cancer_type_det_1
cancer_type_det_2 <- cut_off_dat$cancer_type_det_2
cancer_type_det_3 <- cut_off_dat$cancer_type_det_3
cancer_type_det_4 <- cut_off_dat$cancer_type_det_4


cut_off_dat$cell_1 <- ifelse(cut_off_dat$cell == 'High', 1, 0)
cut_off_dat$cell_2 <- ifelse(cut_off_dat$cell == 'Moderate', 1, 0)

cell_1 <- cut_off_dat$cell_1
cell_2 <- cut_off_dat$cell_2



cut_off_dat$her2_Negative <- ifelse(cut_off_dat$her2 == 'Negative', 1, 0)
#cut_off_dat$her2_Positive <- ifelse(cut_off_dat$her2 == 'Positive', 1, 0)


her2_Negative <- cut_off_dat$her2_Negative
#her2_Positive <- cut_off_dat$her2_Positive



cut_off_dat$Mastectomy <- ifelse(cut_off_dat$type_breast_surg == 'MASTECTOMY', 1, 0)
#cut_off_dat$Breast_Conserving <- ifelse(cut_off_dat$type_breast_surg == 'BREAST CONSERVING', 1, 0)

Mastectomy <- cut_off_dat$Mastectomy
#Breast_Conserving <- cut_off_dat$Breast_Conserving


#cut_off_dat$nhg_1 <- ifelse(cut_off_dat$nhg == 1, 1, 0)
cut_off_dat$nhg_2 <- ifelse(cut_off_dat$nhg == 2, 1, 0)
cut_off_dat$nhg_3 <- ifelse(cut_off_dat$nhg == 3, 1, 0)


#nhg_1 <- cut_off_dat$nhg_1
nhg_2 <- cut_off_dat$nhg_2
nhg_3 <- cut_off_dat$nhg_3


cut_off_dat$Post <- ifelse(cut_off_dat$ims == 'Post', 1, 0)
#cut_off_dat$Pre <- ifelse(cut_off_dat$ims == 'Pre', 1, 0)


Post <- cut_off_dat$Post
#Pre <- cut_off_dat$Pre


#cut_off_dat$gcs_1 <- ifelse(cut_off_dat$gcs == 'ER-/HER2-', 1, 0)
cut_off_dat$gcs_2 <- ifelse(cut_off_dat$gcs == 'ER+/HER2- High Prolif', 1, 0)
cut_off_dat$gcs_3 <- ifelse(cut_off_dat$gcs == 'ER+/HER2- Low Prolif', 1, 0)
cut_off_dat$gcs_4 <- ifelse(cut_off_dat$gcs == 'HER2+', 1, 0)


#gcs_1 <- cut_off_dat$gcs_1
gcs_2 <- cut_off_dat$gcs_2
gcs_3 <- cut_off_dat$gcs_3
gcs_4 <- cut_off_dat$gcs_4


cut_off_dat$cohort_1 <- ifelse(cut_off_dat$cohort == 1, 1, 0)
cut_off_dat$cohort_2 <- ifelse(cut_off_dat$cohort == 2, 1, 0)
cut_off_dat$cohort_3 <- ifelse(cut_off_dat$cohort == 3, 1, 0)
#cut_off_dat$cohort_5 <- ifelse(cut_off_dat$cohort == 5, 1, 0)


cohort_1 <- cut_off_dat$cohort_1
cohort_2 <- cut_off_dat$cohort_2
cohort_3 <- cut_off_dat$cohort_3
#cohort_5 <- cut_off_dat$cohort_5


#cut_off_dat$chemo_yes <- ifelse(cut_off_dat$chemo == 'YES', 1, 0)
cut_off_dat$chemo_no <- ifelse(cut_off_dat$chemo == 'NO', 1, 0)

#chemo_yes <- cut_off_dat$chemo_yes
chemo_no <- cut_off_dat$chemo_no


cut_off_dat$hormone_yes <- ifelse(cut_off_dat$hormone == 'YES', 1, 0)
#cut_off_dat$hormone_no <- ifelse(cut_off_dat$hormone == 'NO', 1, 0)

hormone_yes <- cut_off_dat$hormone_yes
#hormone_no <- cut_off_dat$hormone_no


cut_off_dat$radio_yes <- ifelse(cut_off_dat$radio == 'YES', 1, 0)
#cut_off_dat$radio_no <- ifelse(cut_off_dat$radio == 'NO', 1, 0)

radio_yes <- cut_off_dat$radio_yes
#radio_no <- cut_off_dat$radio_no


cut_off_dat$ptl_left <- ifelse(cut_off_dat$ptl == 'Left', 1, 0)

ptl_left <- cut_off_dat$ptl_left



cut_off_dat$prs_Negative <- ifelse(cut_off_dat$prs == 'Negative', 1, 0)
prs_Negative <- cut_off_dat$prs_Negative



#cut_off_dat$tumor_stage_1 <- ifelse(cut_off_dat$tumor_stage == 1, 1, 0)
cut_off_dat$tumor_stage_2 <- ifelse(cut_off_dat$tumor_stage == 2, 1, 0)
cut_off_dat$tumor_stage_3 <- ifelse(cut_off_dat$tumor_stage == 3, 1, 0)
cut_off_dat$tumor_stage_4 <- ifelse(cut_off_dat$tumor_stage == 4, 1, 0)


#tumor_stage_1 <- cut_off_dat$tumor_stage_1 
tumor_stage_2 <- cut_off_dat$tumor_stage_2
tumor_stage_3 <- cut_off_dat$tumor_stage_3
tumor_stage_4 <- cut_off_dat$tumor_stage_4





##############################################      Creating a data frame with the dummies


dummy_predictors <- data.frame( age = age, Mastectomy = Mastectomy, cancer_type_det_4 = cancer_type_det_4, 
cancer_type_det_2 = cancer_type_det_2, cancer_type_det_3 = cancer_type_det_3, tumor_size = tumor_size, lnep = lnep, 
cell_1 = cell_1, cell_2 = cell_2, chemo_no = chemo_no, cohort_1 = cohort_1, cohort_2 = cohort_2, cohort_3 = cohort_3, nhg_2 = nhg_2, 
nhg_3 = nhg_3, her2_Negative = her2_Negative, hormone_yes = hormone_yes, Post = Post, ptl_left = ptl_left, mc = mc, 
npi = npi, prs_Negative = prs_Negative, radio_yes = radio_yes, gcs_4 = gcs_4, gcs_2 = gcs_2, gcs_3 = gcs_3, Pam50_1 = Pam50_1,
Pam50_2 = Pam50_2,Pam50_3 = Pam50_3, Pam50_4 = Pam50_4, Pam50_5 = Pam50_5, Pam50_6 = Pam50_6, tumor_stage_4 = tumor_stage_4,
tumor_stage_2 = tumor_stage_2, tumor_stage_3 = tumor_stage_3, Time = Time, event = event )


round(apply(dummy_predictors,2,var),2)


##############################################      Penalized regression




library(glmnet)


X <- as.matrix(dummy_predictors[, -c(36,37)])

Y <- matrix(c(Time,event),byrow = F, ncol = 2)

colnames(Y) <- c('time', 'status')






# Standardization should be applied!!!







# We apply the glmnet function to compute the solution path under default settings:
fit <- glmnet(X, Y, family = "cox")




# Each curve corresponds to a variable. It shows the path of its coefficient against the L1-norm of the whole coefficient vector as λ varies.
# The axis above indicates the number of nonzero coefficients at the current λ, which is the effective degrees of freedom (df) for the lasso.
# Users may also wish to annotate the curves: this can be done by setting label = TRUE in the plot command.
plot(fit)


# We can plot the coefficients with the plot method:
plot(fit, label = T)







print(fit)
# It shows from left to right the number of nonzero coefficients (Df), the percent (of null) deviance explained (%dev) and the value of λ
# (Lambda). Although glmnet fits the model for 100 values of lambda by default, it stops early if %dev does not change sufficently from one lambda to the next (typically near the end of the path.) Here we have truncated the prinout for brevity.




set.seed(1)
# The function glmnet returns a sequence of models for the users to choose from. In many cases, users may prefer the software to select one of them.
# Cross-validation is perhaps the simplest and most widely used method for that task. cv.glmnet is the main function to do cross-validation here,
# along with various supporting methods such as plotting and prediction.
cvfit <- cv.glmnet(X, Y, family = "cox", type.measure = "C")



# Once fit, we can view the optimal λ value and a cross validated error plot to help evaluate our model.
# cv.glmnet returns a cv.glmnet object, a list with all the ingredients of the cross-validated fit. As with glmnet, we do not encourage users to 
# extract the components directly except for viewing the selected values of λ.
# The package provides well-designed functions for potential tasks. For example, we can plot the object:
plot(cvfit)



# As with other families, the left vertical line in our plot shows us where the CV-error curve hits its minimum.
# The right vertical line shows us the most regularized model with CV-error within 1 standard deviation of the
# minimum. We also extract such optimal λ’s:
cvfit$lambda.min

cvfit$lambda.1se



coef(cvfit, s = "lambda.min")

coef(cvfit, s = "lambda.lse")













##############################################      Principal Component Analysis

library(survival)
library(survminer)
library(FactoMineR)
library(factoextra)

str(dummy_predictors)


# Standardize the data
#dummy_predictors_std <- dummy_predictors[,-c(36,37)]
#dummy_predictors_std <- scale(dummy_predictors_std)
#dummy_predictors_std <- as.data.frame(dummy_predictors_std)
#str(dummy_predictors_std)


library(stats)

# PCA is an unsupervised, geometric method. Serves as a first step data analysis, which involves multivariate data.
# Can help data exploration by:
# Identifying outliers 
# Possible groups
# Trends
# Collinearities
# The aim is to create linear combinations of the original variables so that the linear combinations variables:
# Are uncorrelated
# Contain as much variability as possible
# Lead to a dimensionality reduction
# Assist in quantification of latent characteristics


set.pr <- princomp(scale(dummy_predictors[,-c(36,37)]))
set.pr # small output


summary(set.pr) # more rich output

# The 1st PC explains the 13.2% of the total variability
# The first 2 PCs explain the 21.4% of the total variability
# The first 3 PCs explain the 28.3% of the total variability
# The first 4 PCs explain the 34.3% of the total variability
# We observe that the standard deviations with value more than 1, are the first 14 ones


set.pr$loadings

range(set.pr$scores[,1])
range(set.pr$scores[,2])

screeplot(set.pr,type = "lines")
# The screeplot indicates to use the first 2 PCs as the greatest angle is located between the first and the second 
# component.


# Biplot (or Collinearity plot)
biplot(set.pr, choices = c(1,2), scale = 1, cex = 0.6, pc.biplot = F, col = c('blue','red'), xlab = 'First Component', ylab = 'Second Component'
       )
# Red color: The loadings
# position of each variable + loading -> inference

# Positive loadings indicate that a variable and a principal component are positively correlated whereas negative loadings
# indicate a negative correlation. When loadings are large, (+ or -) it indicates that a variable has a strong effect on
# that principal component.


# Possible outliers when looking the biplot from PC1,PC2
cut_off_dat[172,] # observation 172 has tumor_size 180 while the mean tumor_size is 26 (and reasonably has the MAX npi!)
cut_off_dat[318,] # observation 318 has tumor_size 160 while the mean tumor_size is 26 (and reasonably has VERY high npi!)
cut_off_dat[900,] # observation 900 has npi 6.16 while the mean npi is 4.129, lnpe 16 while the mean lnep is 1.889, tumor_size 80 while the mean tumor_size is 26
cut_off_dat[920,] # observation 920 has tumor_size 2 while the mean tumor_size is 26
cut_off_dat[380,] # observation 380 has tmb in Q1, tumor_size less than Q1, npi less than Q1, lnep 0 (while the mean is 1.889 )
cut_off_dat[527,] # observation 527 has the min age (21.93), tumor_size 35, lnep 6 (while the mean is 1.889), tumor_stage 3(?!)
cut_off_dat[444,] # observation 444 has tumor_size less than Q1, lnep 1 (while the mean is 1.889), mc less than Q1, tmb less than Q1
cut_off_dat[954,] # observation 954 has tumor_size 60 while the mean is 26, lnep 4 (while the mean is 1.889), npi 6.12 (near the max which is 6.360)
cut_off_dat[135,] # observation 135 has the max lnep, 41 while the mean lnep is 1.889, npi close to the max npi, 
cut_off_dat[159,] # observation 159 has tumor_size 50 (while the mean tumor_size is 26), lnep 19 (while the mean lnep is 1.889), the min mc, the min tmb, npi 6.1 (close to the max)



# Possible outliers when looking the biplot from PC1,PC3
biplot(set.pr, choices = c(1,3), scale = 1, cex = 0.7, pc.biplot = F , col = c('blue','red'), xlab = 'First Component', ylab = 'Third Component')

cut_off_dat[247,] # observation 172 has tumor_size 65 (while the mean tumor_size is 26), npi 6.13 (mean npi is 4.1), lnep 23 (mean lnep is 1.8)
cut_off_dat[677,] # observation 677 has the max mc (and the max tmb), the min lnep
cut_off_dat[135,] # observation 135 has the max lnep (and quite high npi)
cut_off_dat[318,] # observation 318 has tumor_size 160, npi 6.32, lnep 22



# Possible outliers when looking the biplot from PC2,PC3
biplot(set.pr, choices = c(2,3), scale = 1, cex = 0.7, pc.biplot = F , col = c('blue','red'), xlab = 'Second Component', ylab = 'Third Component')

cut_off_dat[247,] # observation 172 has tumor_size 65 (while the mean tumor_size is 26), npi 6.13 (mean npi is 4.1), lnep 23 (mean lnep is 1.8)
cut_off_dat[172,] # observation 172 has tumor_size 180 while the mean tumor_size is 26 (and reasonably has the MAX npi!)
cut_off_dat[677,] # observation 677 has the max mc (and the max tmb), the min lnep, 
cut_off_dat[615,] # observation 615 has 41 mc (the mean mc is 5.4) and close to the max tmb
cut_off_dat[14,] # observation 14 has 6.3 npi (max npi is 6.36), tumor_size 150 (while the mean tumor_size is 26)



# Possible outliers when looking the biplot from PC1,PC4
biplot(set.pr, choices = c(1,4), scale = 1, cex = 0.7, pc.biplot = F, col = c('blue','red'), xlab = 'First Component', ylab = 'Fourth Component' )

cut_off_dat[900,] # observation 900 has npi 6.16 while the mean npi is 4.129, lnpe 16 while the mean lnep is 1.889, tumor_size 80 while the mean tumor_size is 26
cut_off_dat[969,] # observation 969 has age 90.08, tumor stage 3
cut_off_dat[135,] # observation 135 has the max lnep (and quite high npi)
cut_off_dat[954,] # observation 954 has tumor_size 60 while the mean is 26, lnep 4 (while the mean is 1.889), npi 6.12 (near the max which is 6.360)
cut_off_dat[318,] # observation 318 has tumor_size 160, npi 6.32, lnep 22




# Possible outliers when looking the biplot from PC2,PC4
biplot(set.pr, choices = c(2,4), scale = 1, cex = 0.7, pc.biplot = F )

cut_off_dat[969,] # observation 969 has age 90.08, tumor stage 3
cut_off_dat[900,] # observation 900 has npi 6.16 while the mean npi is 4.129, lnpe 16 while the mean lnep is 1.889, tumor_size 80 while the mean tumor_size is 26
cut_off_dat[14,] # observation 14 has 6.3 npi (max npi is 6.36), tumor_size 150 (while the mean tumor_size is 26)



# Possible outliers when looking the biplot from PC3,PC4
biplot(set.pr, choices = c(3,4), scale = 1, cex = 0.7, pc.biplot = F )

cut_off_dat[969,] # observation 969 has age 90.08, tumor stage 3
cut_off_dat[900,] # observation 900 has npi 6.16 while the mean npi is 4.129, lnpe 16 while the mean lnep is 1.889, tumor_size 80 while the mean tumor_size is 26
cut_off_dat[247,] # observation 172 has tumor_size 65 (while the mean tumor_size is 26), npi 6.13 (mean npi is 4.1), lnep 23 (mean lnep is 1.8)
cut_off_dat[14,] # observation 14 has 6.3 npi (max npi is 6.36), tumor_size 150 (while the mean tumor_size is 26)
cut_off_dat[677,] # observation 677 has the max mc (and the max tmb), the min lnep, 
cut_off_dat[615,] # observation 615 has 41 mc (the mean mc is 5.4) and close to the max tmb


## Common (2+) observations:
##
# 172 (2 times), 318 (3 times), 247 (3 times), 135 (3 times), 900 (4 times), 677 (2 times), 954 (2 times), 969 (3 times), 14 (3 times), 615 (2 times)
# rule to exclude: 2+ or 3+







# Graph of variables. Positive correlated variables point to the same side of the plot. Negative correlated variables
# point to opposite sides of the graph.
#fviz_pca_var(set.pr,
#             col.var = "contrib", # Color by contributions to the PC
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE # Avoid text overlapping
#)






####################################           Final predictors matrix

age_star <- c("Age")
cancer_type_det_star <- c("Cancer Type Detailed")
lnep_star <- c("Lymph Nodes Examined Positive")
gcs_star <- c("3 Gene Classifier Subtype")
tumor_size_star <- c("Tumor Size")


var_sel_mat <- matrix( c(age_star,
                         cancer_type_det_star,
                        lnep_star,
                        gcs_star, tumor_size_star) , 
                       nrow = 5, ncol = 1
)

var_sel_mat




colnames(var_sel_mat) <- c("METABRIC breast cancer predictors")
rownames(var_sel_mat) <- c("1","2","3","4","5")
var_sel_mat




library(kableExtra)

# Convert matrix to a data frame
var_sel_df <- as.data.frame(var_sel_mat)




# Create a fancy matrix with the mean survival time estimates for the 50% scenario
caption.2.1 <- "<b>Final predictors selected for the modeling</b>"
# Create a styled table using kableExtra
table.2.1 <- kable(var_sel_df, format = "html", table.attr = "class='table'", 
                   caption = caption.2.1) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Display the table
table.2.1





####################################           Modeling via survHE


# formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + nhg # Runs via MLE only! Doesn't run via HMC
# Note that gcs and nhg are NOT independent!

# formula <- Surv(Time, event) ~ age + tumor_size + lnep + cancer_type_det + tumor_stage + npi # Does NOT run
# formula <- Surv(Time, event) ~ age + tumor_size + lnep + cancer_type_det + tumor_stage + her2 + cohort # Does NOT run

# formula <- Surv(Time, event) ~ age + tumor_size + lnep + cancer_type_det + tumor_stage + her2 + ims # Does NOT run

# formula <- Surv(Time, event) ~ age + tumor_size + lnep + cancer_type_det + tumor_stage + her2 # Runs


formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 +
cancer_type_det_4 ##### Runs!

formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + tumor_stage ### Runs

formula2 <- Surv(Time, event) ~ age + tumor_size + lnep + gcs_1 + gcs_2 + gcs_3 + cancer_type_det_1 + cancer_type_det_2 +
+ cancer_type_det_3 + tumor_stage_1 + tumor_stage_2 + tumor_stage_3

# formula <- Surv(Time, event) ~ age  + cancer_type_det + tumor_size + lnep + gcs + tumor_stage + ims # Does NOT run

# formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + her2 # Runs
# Note that gcs and her2 are NOT independent

# formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det # Runs


# formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + npi # Does NOT run
# formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + her2 + ims # Does NOT run
# formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + her2 + cohort # Does NOT run


# formula <- Surv(Time, event) ~ age + tumor_size + lnep + cancer_type_det + npi + her2  # Does NOT run
# formula <- Surv(Time, event) ~ age + tumor_size + lnep + cancer_type_det + npi # Does NOT run


# Older versions of the formula i tried:
# formula <- Surv(Time, event) ~ age + type_breast_surg + cohort + nhg + ims + npi + gcs # runs
# formula <- Surv(Time, event) ~ age + type_breast_surg + cohort + nhg + ims + npi + gcs + radio # runs
# formula <- Surv(Time, event) ~ age + type_breast_surg + cohort + nhg + ims + npi + gcs + radio + hormone + chemo  # Does NOT run
# formula <- Surv(Time, event) ~ age + type_breast_surg + cohort + nhg + ims + npi + gcs + radio + hormone #  Does not run!
# formula <- Surv(Time, event) ~ age + type_breast_surg + cohort + nhg + ims + npi + gcs + radio + chemo  # Does not run!


library(survHE)
## run_mle
# Run a selection of models using the fit.models function in survHE
# First, defines the vector of models to be used
mods <- c("exponential","weibull","loglogistic",'lognormal','gompertz','gamma','gengamma')


gc()

# Then the formula specifying the linear predictor:
# formula = Surv(Time, event) ~ age + cancer_type_det + tumor_size + nhg + her2 + hormone + ims + lnep + npi + radio + gcs 
# 11 from the automatic procedures



# Run the models using MLE via flexsurv & show the elements of the resulting survHE object
m1 <- fit.models(formula = formula, data = cut_off_dat, distr = mods)


gc()
 



# Explore the output
class(m1)
names(m1)
names(m1$models)
names(m1$models[[1]])
names(m1$models$Exponential)
m1$model.fitting
m1$method



## Visualising the results
# Use the print method for survHE objects. This uses the first model stored in m1, by default
print(m1)

# But can choose to show the fifth model
print(m1, mod = 5)

# Or use the original flexsurv table instead
print(m1, mod = 1, original = TRUE, digits = 3)



# Open a new window to plot the next graph
windows()

# Another important issue is to do model fitting - we can use survHE 
model.fit.plot(m1)

# Selecting the type of statistic ("aic" vs "bic" vs "dic")
model.fit.plot(m1, type = "bic")


# Adjust font size globally for all plots
par(cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)



## run_hmc
# NB: it is possible to use the option 'refresh = 0' to prevent rstan from printing the progress indicator

formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + tumor_stage # AIC = 6221.229

formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + nhg # Does NOT run
formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det # AIC = 6223.740
formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs # AIC = 6223.209

mods <- c("exponential","weibull","loglogistic","lognormal","gompertz")




m3 <- fit.models(formula, data = cut_off_dat, distr = mods, "hmc")

print(m3) 


## Traceplots
rstan::traceplot(m3$models[[2]])


## ACF
rstan::stan_ac(m3$models[[2]])






## rps
# We can also do "advanced" models, i.e. Royston-Parmar splines & Poly-Weibull
# RPS is available under MLE & Bayesian HMC modelling:

formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + tumor_stage 

# old formula = Surv(Time,event) ~ age + cancer_type_det + tumor_size + lnep + gcs + ims

m6 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 2)
print(m6) # AIC = 6165.091/ With the older formula: AIC = 6164.651


mod7 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 3)
print(mod7) # AIC = 6162.475/ With the older formula: AIC =  6161.971


mod8 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 4)
print(mod8) # AIC = 6152.552/ With the older formula: AIC = 6151.890


mod9 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 5)
print(mod9) # AIC = 6155.458 / With the older formula: AIC = 6154.809


mod10 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 6)
print(mod10) # AIC = / With the older formula: AIC = 6154.234


mod11 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 7)
print(mod11) # AIC = / With the older formula: AIC = 6156.710


mod12 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 8)
print(mod12) # AIC = / With the older formula: AIC = 6156.162


mod13 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 9)
print(mod13) # AIC = / With the older formula: AIC = 6158.623


mod14 <- fit.models(formula = formula, data = cut_off_dat, distr = "rps", k = 10)
print(mod14) # AIC = / With the older formula:  AIC = 6163.178









## Poly-weibull

# For the Poly-Weibull, we need to specify a list of formulae before running the 
# Bayesian/HMC model
# k = 2

# with formula.pw <- list(formula <- Surv(Time, event) ~ age + cancer_type_det + tumor_size + lnep + gcs + ims, 
                        # formula <- Surv(Time, event) ~ age + cancer_type_det + tumor_size + lnep + gcs + ims) 
# AIC = 6198.707

formula.pw <- list(formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + tumor_stage,
                   formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + tumor_stage)
mod15  <- poly.weibull(formula.pw, cut_off_dat, cores = 4)


print(mod15) # AIC = 6202.800



formula.pw <- list(formula <- Surv(Time, event) ~ age + tumor_size + lnep + cancer_type_det + tumor_stage + her2,
                   formula <- Surv(Time, event) ~ age + tumor_size + lnep + cancer_type_det + tumor_stage + her2)
mod15  <- poly.weibull(formula.pw, cut_off_dat, cores = 4)




# k = 3
# with formula.pw <- list(formula <- Surv(Time, event) ~ age + cancer_type_det + tumor_size + lnep + gcs + ims , 
                   # formula <- Surv(Time, event) ~ age + cancer_type_det + tumor_size + lnep + gcs + ims ,
                   # formula <- Surv(Time, event) ~ age + cancer_type_det + tumor_size + lnep + gcs + ims)
# AIC = 6222.382

formula.pw <- list(formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + tumor_stage,
                   formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + tumor_stage,
                   formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs + cancer_type_det + tumor_stage)

mod16  <- poly.weibull(formula.pw, cut_off_dat, cores = 4)


print(mod16) # AIC = 6230.730









##############################################      Mean Survival Time Confirmation


library(cubature)

# Define the integrand function 
integrand <- function(x) {
  S <- exp(-(lambda1 * x)) # The Exponential Survivor Function
  return(S)
}

# Evaluate the integral
result <- adaptIntegrate(integrand, lowerLimit = 0, upperLimit = Inf)$integral
result







##############################################     Plotting the KM of the cut_off_dat up to 400 months




# KM for the cut off data
km <- survfit(Surv(Time, event) ~ 1, type = "kaplan-meier", data = data)

# Plot the KM while expanding the xlim to 400 months
plot(km, conf.int = F, lwd = 3, xlab = "Months", main = 'Kaplan-Meier plot', ylab = 'S(t)', xlim = c(0,400), las = 1 ) 



library(ggfortify)
autoplot(km, lwd = 2, xlab = "Months", main = 'Kaplan-Meier plot of the full data', ylab = 'S(t)', las = 1 )
abline(h = 0.125, col = 'darkred',lty =2)



autoplot(km.full, lwd = 2, xlab = "Months", main = 'Kaplan-Meier Plot of the full data', ylab = 'S(t)', xlim = c(0,400), las = 1 )






##############################################      Extrapolation via survival package 


formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 + 
  cancer_type_det_4


# Fit an Exponential model via survival package
model.expo <- survreg(Surv(Time, event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + 
                        cancer_type_det_3 + cancer_type_det_4, data = cut_off_dat, dist = "exp")
print(model.expo)

# For the Exponential distribution: Assuming that h(t) = λ = exp(−Intercept), we calculate λ by:
lambda1 <- exp(-(model.expo$coefficients[1] + model.expo$coefficients[2]*mean(age) + model.expo$coefficients[3]*mean(tumor_size) +
model.expo$coefficients[4]*mean(lnep) + model.expo$coefficients[5]*mean(gcs_2) + model.expo$coefficients[6]*mean(gcs_3) +
model.expo$coefficients[7]*mean(gcs_4) + model.expo$coefficients[8]*mean(cancer_type_det_2) + 
model.expo$coefficients[9]*mean(cancer_type_det_3) + model.expo$coefficients[10]*mean(cancer_type_det_4 )))


t <- seq(161,400,1)
St2 <- exp(-(lambda1 * t))
ex <- as.data.frame(cbind(t = t, St = St2))
summary(ex$St)


# Adding it to the plot
lines(ex$t, ex$St, lty = 1, col = "blue", lwd = 3)








# Second way to build the model and plot the extrapolated curve (via survHE)
print(m1$models$Exponential) # Exponentiated coeffs (not needed)


t <- seq(161,400,1) # Create a vector of time t


m1$models$Exponential$coefficients # Extract the (log) coeffs


# Calculate the rate coefficient
lambda2 <- exp( m1$models$Exponential$coefficients[1] + m1$models$Exponential$coefficients[2]*mean(age) +
                  m1$models$Exponential$coefficients[3]*mean(tumor_size) + m1$models$Exponential$coefficients[4]*mean(lnep) +
                  m1$models$Exponential$coefficients[5]*mean(gcs_2) + m1$models$Exponential$coefficients[6]*mean(gcs_3) +
                  m1$models$Exponential$coefficients[7]*mean(gcs_4) + m1$models$Exponential$coefficients[8]*mean(cancer_type_det_2) +
                  m1$models$Exponential$coefficients[9]*mean(cancer_type_det_3) + m1$models$Exponential$coefficients[10]*mean(cancer_type_det_4)
                  )

# Create a vector to store the values of the survivor function
expo_survivor_func <- numeric(length(t)) 

# Survivor function for the exponential model
for (i in 1:length(t)){
  expo_survivor_func[i] <- exp(1)^(-lambda2*t[i])
}
summary(expo_survivor_func)


# Adding it to the plot
lines(t, expo_survivor_func, lty = 1, col = "pink", lwd = 2)









## 3rd way to build the survivor function for the exponential parametric survival model

# Cumulative hazard function for the exponential distribution
# cum_haz_expo_func <- numeric(length(t))
# for (i in 1:length(t)) {
#   cum_haz_expo_func[i] <- lambda*t[i]
# }
#summary(cum_haz_expo_func)


# Survivor function for the exponential model (via the expression that connects the cumulative hazard and the survivor function)
#surv_expo <- numeric(length(t))
#for (i in 1:length(t)) {
#  surv_expo[i] <- exp(-cum_haz_expo_func[i])
#}
#summary(surv_expo) 






##################################            Extrapolation via rms package



library(rms) # Regression Modeling Strategies

# psm is a modification of Therneau's survreg function for fitting the accelerated failure time family of parametric 
# survival models. psm uses the rms class for automatic anova, fastbw, calibrate, validate, and other functions. 
# Hazard.psm, Survival.psm, Quantile.psm, and Mean.psm create S functions that evaluate the hazard, survival, quantile,
# and mean (expected value) functions analytically, as functions of time or probabilities and the linear predictor 
# values. The Nagelkerke R^2 and and adjusted Maddala-Cox-Snell R^2 are computed. For the latter the notation is 
# R2(p,m) where p is the number of regression coefficients being adjusted for and m is the effective sample size 
# (number of uncensored observations). See R2Measures for more information.

d <- datadist(cut_off_dat)
options(datadist = "d")

formula <- Surv(Time, event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 + cancer_type_det_4


# Fit an Exponential model via the rms package
expo <- psm((Surv(Time,event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 + cancer_type_det_4 ), 
            dist = 'exponential', data = cut_off_dat)
print(expo, title = 'Exponential Model')


# Computes predicted survival probabilities or hazards and optionally confidence limits (for survival only) for 
# parametric survival models fitted with psm. If getting predictions for more than one observation, times must be 
# specified. For a model without predictors, no input data are specified.
probs_expo <- survest(expo, cut_off_dat, times = seq(161,400,1), conf.int = FALSE, what = 'survival')
head(probs_expo)


mean_surv_month_expo <- apply(probs_expo,2,mean)

summary(mean_surv_month_expo)

t <- seq(161,400,1)

lines(t,mean_surv_month_expo,col='red',lwd = 2)





# Fit a Weibull model via the rms package
wei.rms <- psm((Surv(Time,event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 + cancer_type_det_4 ), 
           dist = 'weibull', data = cut_off_dat)


print(wei.rms, title = 'Weibull Model')

probs_wei.rms <- survest(wei.rms, cut_off_dat, times = seq(161,400,1), conf.int = FALSE, what = 'survival')
head(probs_wei.rms)


mean_surv_month_wei.rms <- apply(probs_wei.rms,2,mean)
summary(mean_surv_month_wei.rms)

lines(t, mean_surv_month_wei.rms, lty = 1, col = 'orange', lwd = 2)







####################################         Extrapolation via survival package


# Fit a Weibull AFT model via survival package
wei_ <- survreg( Surv(Time,event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 +
                   cancer_type_det_4,
                 data = cut_off_dat, dist = "weibull")
print(wei_)

# Using the Weibull distribution mu (the intercept) = -log(lambda) and sigma (scale) = 1/p
# Thus, the scale parameter: lambda = exp(-mu) and p = 1/sigma
p11 <- 1/wei_$scale
p11

# Calculate the shape
lambda33 <- exp(-(wei_$coefficients[1] + wei_$coefficients[2]*mean(age) + wei_$coefficients[3]*mean(tumor_size)
                + wei_$coefficients[4]*mean(lnep) + wei_$coefficients[5]*mean(gcs_2) + wei_$coefficients[6]*mean(gcs_3)
                + wei_$coefficients[7]*mean(gcs_4) + wei_$coefficients[8]*mean(cancer_type_det_2)
                + wei_$coefficients[9]*mean(cancer_type_det_3) + wei_$coefficients[10]*mean(cancer_type_det_4) ))

lambda33


St <- exp(-(lambda33*t)^p11)
weib <- as.data.frame(cbind(t=t,St=St))
head(weib)
summary(weib$St)


lines(weib$t, weib$St, lty = 1, col = 'purple', lwd = 3)



# Create a legend
legend('topright', legend=c("Exponential - Survival/survHE", "Exponential - RMS", "Weibull (AFT) - RMS",'Weibull (AFT) - Survival/survHE'),
       col = c("blue","red","orange","purple"), lwd = 2, cex = 0.8)






# Third way to build the Weibull AFT model and plot the extrapolated curve via survHE
print(m1$models$`Weibull (AFT)`) # Exponentiated coeffs (not needed)


m1$models$`Weibull (AFT)`$coefficients # Extract the (log) coeffs


pp <- exp(m1$models$`Weibull (AFT)`$coefficients[1]) # Calculate the shape
pp


# Calculate the scale coefficient
lambda4 <- exp( -(m1$models$`Weibull (AFT)`$coefficients[2] + m1$models$`Weibull (AFT)`$coefficients[3]*mean(age) +
                  m1$models$`Weibull (AFT)`$coefficients[4]*mean(tumor_size) + m1$models$`Weibull (AFT)`$coefficients[5]*mean(lnep) +
                  m1$models$`Weibull (AFT)`$coefficients[6]*mean(gcs_2) + m1$models$`Weibull (AFT)`$coefficients[7]*mean(gcs_3) +
                  m1$models$`Weibull (AFT)`$coefficients[8]*mean(gcs_4) + m1$models$`Weibull (AFT)`$coefficients[9]*mean(cancer_type_det_2) +
                  m1$models$`Weibull (AFT)`$coefficients[10]*mean(cancer_type_det_3) + m1$models$`Weibull (AFT)`$coefficients[11]*mean(cancer_type_det_4)
))
lambda4

# Create a vector to store the values of the survivor function
wei_survivor_func <- numeric(length(t))

# Survivor function for the exponential model
for (i in 1:length(t)){
  wei_survivor_func[i] <- exp(-(lambda4*t[i])^pp )
}
summary(wei_survivor_func) # Exactly the same as the Weibull AFT build with survival package



# Compare the summaries

# Expo with 5 predictors (both numerical and categorical)
summary(ex$St) # survreg
summary(expo_survivor_func) # survHE
summary(mean_surv_month_expo) # Rms


# Weibull (AFT) with 5 predictors (both numerical and categorical)
summary(weib$St) # survreg 
summary(mean_surv_month_wei.rms) # rms
summary(wei_survivor_func) # survHE






####################################         Fitting  Gompertz model via flexsurv package 


model.gomp <- flexsurvreg(Surv(Time, event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 +
                            cancer_type_det_4, data = cut_off_dat, dist = "gompertz")
print(model.gomp)

print(m1$models$Gompertz) # exponentiated coeffs

m1$models$Gompertz$coefficients # Extract the (log) coeffs



# rate parameter, MUST BE POSITIVE!
b2 <- exp( m1$models$Gompertz$coefficients[2] + m1$models$Gompertz$coefficients[3]*mean(age) +
       m1$models$Gompertz$coefficients[4]*mean(tumor_size) + m1$models$Gompertz$coefficients[5]*mean(lnep) +
       m1$models$Gompertz$coefficients[6]*mean(gcs_2) + m1$models$Gompertz$coefficients[7]*mean(gcs_3) +
       m1$models$Gompertz$coefficients[8]*mean(gcs_4) + m1$models$Gompertz$coefficients[9]*mean(cancer_type_det_2) +
       m1$models$Gompertz$coefficients[10]*mean(cancer_type_det_3) + m1$models$Gompertz$coefficients[11]*mean(cancer_type_det_4 )     ) # Calculate the rate parameter
b2

# shape parameter 
a2 <- m1$models$Gompertz$coefficients[1]
a2
# Note that a < 0 is permitted, in which case S(t) tends to a non-zero probability as t increases, i.e. a probability of living forever.


# Create a vector to store the values of the survivor function
gomp_survivor_func <- numeric(length(t)) 


# Survivor function for the Gompertz distribution
for (i in 1:length(t)){
  gomp_survivor_func[i] <- exp( -(b2/a2)*(exp(a2*t[i])-1))
}
summary(gomp_survivor_func)


lines(t, gomp_survivor_func, lty = 1, col='cyan', lwd = 3)

# Create a legend
legend('topright', legend=c("Exponential - Survival/survHE", "Exponential - RMS", "Weibull (AFT) - RMS",
                            'Weibull (AFT) - Survival/survHE',"Gompertz - Survival/survHE"),
       col = c("blue","red","orange","purple","cyan"), lwd = 2, cex = 1)







####################################         Fitting  Log-Normal model via flexsurv package


model.log.norm <- flexsurvreg(Surv(Time, event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 +
                            cancer_type_det_4, data = cut_off_dat, dist = "lognorm")
print(model.log.norm)

#plot(model.log.norm, type = 'survival')


print(m1$models$`log-Normal`) # exponentiated coeffs

m1$models$`log-Normal`$coefficients # Extract the (log) coeffs


log_norm <- survreg(Surv(Time, event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 +
                      cancer_type_det_4, data = cut_off_dat, dist = "lognorm")

summary(log_norm)

# shape parameter
s <- 1/log_norm$scale
s


# scale parameter
mu <- exp( -(m1$models$`log-Normal`$coefficients[1] + m1$models$`log-Normal`$coefficients[3]*mean(age) +
             m1$models$`log-Normal`$coefficients[4]*mean(tumor_size) + m1$models$`log-Normal`$coefficients[5]*mean(lnep) +
             m1$models$`log-Normal`$coefficients[6]*mean(gcs_2) + m1$models$`log-Normal`$coefficients[7]*mean(gcs_3) +
             m1$models$`log-Normal`$coefficients[8]*mean(gcs_4) + m1$models$`log-Normal`$coefficients[9]*mean(cancer_type_det_2) +
             m1$models$`log-Normal`$coefficients[10]*mean(cancer_type_det_3) + m1$models$`log-Normal`$coefficients[11]*mean(cancer_type_det_4 )  )   ) # Calculate the rate parameter
mu






# Survivor function for the log normal model via the survHE package
x <- seq(161,400,1)
S.log.norm <- 1-pnorm(log(x), log_norm$coef[1] + log_norm$coef[2]*mean(age) + log_norm$coef[3]*mean(tumor_size) + log_norm$coef[4]*mean(lnep) +
                        log_norm$coef[5]*mean(gcs_2) + log_norm$coef[6]*mean(gcs_3) + log_norm$coef[7]*mean(gcs_4) +  
                        log_norm$coef[8]*mean(cancer_type_det_2) + log_norm$coef[9]*mean(cancer_type_det_3) + log_norm$coef[10]*mean(cancer_type_det_4) 
                      
                      , log_norm$scale)



lines(t, S.log.norm, lty = 1, col = 'red', lwd = 3)



legend('topright', legend=c("Exponential", "Weibull (AFT)",
                            "Gompertz", "Log-Normal",
                            "Log-Logistic"),
       col=c("blue","purple","cyan","red","green"), lwd = 2, cex = 1.1)


abline(h = 0.5, lty = 2, col = "darkred")







####################################         Fitting a Log-Logistic model via the survival package


log.logistic <- survreg( Surv(Time, event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + 
            gcs_4 + cancer_type_det_2 + cancer_type_det_3 + cancer_type_det_4 , data = cut_off_dat, dist = "loglogistic")
print(log.logistic)



# Using the Log-logistic distribution mu (the intercept) = -log(lambda) and sigma (scale) = 1/p
# Thus, the scale parameter: lambda = exp(-mu) and p = 1/sigma
scale <- 1/log.logistic$scale

lambda5 <- exp(-(log.logistic$coefficients[1] + log.logistic$coefficients[2]*mean(age) + log.logistic$coefficients[3]*mean(tumor_size)
+ log.logistic$coefficients[4]*mean(lnep) + log.logistic$coefficients[5]*mean(gcs_2) + log.logistic$coefficients[6]*mean(gcs_3)
+ log.logistic$coefficients[7]*mean(gcs_4) + log.logistic$coefficients[8]*mean(cancer_type_det_2)
+ log.logistic$coefficients[9]*mean(cancer_type_det_3) + log.logistic$coefficients[10]*mean(cancer_type_det_4) ))


St <- 1/( 1 + (lambda5*t)^scale )
log.logis <- as.data.frame(cbind( t = t, St = St ))
head(log.logis)

summary(log.logis$St)


lines(log.logis$t, log.logis$St, lty = 1, col = 'green', lwd = 3)



#legend('topright', legend=c("Exponential - Survival/survHE", "Exponential - RMS", "Weibull (AFT) - RMS",
#                            "Weibull (AFT) - Survival/survHE","Gompertz - Survival/survHE","Log-Logistic - Survival/survHE"),
#       col=c("blue","red","orange","purple","cyan","green"), lwd = 2, cex = 0.8)


#legend('topright', legend=c(
#  "Exponential - Survival/survHE", "Exponential - RMS", "Weibull (AFT) - RMS",
#  "Weibull (AFT) - Survival/survHE", "Gompertz - Survival/survHE", "Log-Logistic - Survival/survHE"
#), col=c("blue","red","orange","purple","cyan","green"), lwd = 2, cex = 0.8, ncol = 2)




#legend('topright', legend=c("Exponential","Weibull (AFT)",
#                          "Gompertz","Log-Normal","Log-Logistic","Kaplan-Meier"),
#       col=c("blue","purple","cyan","red","green","black"), lwd = 3, cex = 0.75, )



#legend('topright', legend = c("Exponential", "Weibull (AFT)", "Gompertz", "Log-Normal", "Log-Logistic", "Kaplan-Meier"),
#       col = c("blue", "purple", "cyan", "red", "green", "black"), lwd = 8, cex = 1.05, ncol = 2)
#

legend('topright', legend = c("Exponential", "Weibull (AFT)", "Gompertz", "Log-Normal", "Log-Logistic", "Kaplan-Meier"),
       col = c("blue", "purple", "cyan", "red", "green", "black"), lwd = 8, cex = 1.7, ncol = 1,
       x.intersp = 0.5, y.intersp = 0.3)









####################################         Fitting a Log-Logistic model via the rms package


log.log.rms <- psm((Surv(Time,event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 + cancer_type_det_4 ), 
               dist = 'loglogistic', data = cut_off_dat)


print(log.log.rms, title = 'Log-Logistic Model')

probs_log.log.rms <- survest(log.log.rms, cut_off_dat, times = seq(161,400,1), conf.int = FALSE, what='survival')
head(probs_log.log.rms)


mean_surv_month_log.log.rms <- apply(probs_log.log.rms,2,mean)
summary(mean_surv_month_log.log.rms)

lines(t, mean_surv_month_log.log.rms, lty = 1, col='lightblue', lwd = 2)



legend('topright', legend=c("Exponential - Survival/survHE", "Exponential - RMS", "Weibull (AFT) - RMS",
                            "Weibull (AFT) - Survival/survHE", "Gompertz - Survival/survHE",
                            "Log-Logistic - Survival/survHE","Log-Logistic - RMS"),
       col=c("blue","red","orange","purple","cyan","green","lightblue"), lwd = 2, cex = 0.8)






####################################         Fitting a Log-Normal model via the rms package


log.norm.rms <- psm((Surv(Time,event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 + cancer_type_det_4 ), 
               dist = 'lognormal', data = cut_off_dat)


print(log.norm.rms, title = 'Log-Normal Model')

probs_log.norm.rms <- survest(log.norm.rms, cut_off_dat, times = seq(161,400,1), conf.int = FALSE, what='survival')
head(probs_log.norm.rms)


mean_surv_month_log.norm.rms <- apply(probs_log.norm.rms,2,mean)
summary(mean_surv_month_log.norm.rms)

lines(t, mean_surv_month_log.norm.rms, lty = 1, col = 'darkred', lwd = 2)


legend('topright', legend = c("Exponential - Survival/survHE", "Exponential - RMS", "Weibull (AFT) - RMS",
                            "Weibull (AFT) - Survival/survHE", "Gompertz - Survival/survHE",
                            "Log-Logistic - Survival/survHE", "Log-Logistic - RMS","Log-Normal - RMS"),
       col=c("blue","red","orange","purple","cyan",
                   "green","lightblue","darkred"), lwd = 2, cex = 0.8)








####################################         Fitting a Normal model via the rms package


gauss.rms <- psm((Surv(Time,event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 + cancer_type_det_4 ), 
                    dist = 'gaussian', data = cut_off_dat)


print(gauss.rms, title = 'Normal Model')

probs_gauss.rms <- survest(gauss.rms, cut_off_dat, times = seq(161,400,1), conf.int = FALSE, what='survival')
head(probs_gauss.rms)


mean_surv_month_gauss.rms <- apply(probs_gauss.rms,2,mean)
summary(mean_surv_month_gauss.rms)

lines(t, mean_surv_month_gauss.rms, lty = 1, col='darkslategrey', lwd=2)


legend('topright', legend=c("Exponential - Survival/survHE", "Exponential - RMS", "Weibull (AFT) - RMS",
                            "Weibull (AFT) - Survival/survHE", "Gompertz - Survival/survHE",
                            "Log-Logistic - Survival/survHE","Log-Logistic - RMS","Log-Normal - RMS","Normal - RMS"),
       col=c("blue","red","orange","purple","cyan","green","lightblue","darkred","darkslategrey"), lwd = 2, cex = 0.8)







####################################         Fit a Logistic model via the rms package


logistic.rms <- psm((Surv(Time,event) ~ age + tumor_size + lnep + gcs_2 + gcs_3 + gcs_4 + cancer_type_det_2 + cancer_type_det_3 + cancer_type_det_4 ), 
                 dist = 'logistic', data = cut_off_dat)


print(logistic.rms, title = 'Logistic Model')

probs_logistic.rms <- survest(logistic.rms, cut_off_dat, times = seq(161,400,1), conf.int = FALSE, what='survival')
head(probs_logistic.rms)


mean_surv_month_logistic.rms <- apply(probs_logistic.rms,2,mean)
summary(mean_surv_month_logistic.rms)

lines(t, mean_surv_month_logistic.rms, lty = 1, col='deeppink', lwd = 2)


legend('topright', legend=c("Exponential - Survival/survHE", "Exponential - RMS", "Weibull (AFT) - RMS",
                            "Weibull (AFT) - Survival/survHE","Gompertz - Survival/survHE",
                            "Log-Logistic - Survival/survHE","Log-Logistic - RMS",
                            "Log-Normal - RMS","Normal - RMS","Logistic - RMS"),
       col=c("blue","red","orange","purple","cyan",
                   "green","lightblue","darkred","darkslategrey","deeppink"), lwd = 2, cex = 0.8)












# Extract the AIC for each model we fitted 
expo_aic <- AIC(model.expo)
expo_aic # AIC = 6221.661

wei_aic <- AIC(wei.rms)
wei_aic # AIC = 6178.101


log_log_aic <- AIC(log.log.rms)
log_log_aic # AIC = 6154.795


log_normal_aic <- AIC(log.norm.rms)
log_normal_aic # AIC = 6176.444

gomp_aic <- AIC(model.gomp)
gomp_aic # AIC = 6204.352



# normal_aic <- extractAIC(gauss.rms)[2]
# normal_aic # AIC = 6412.939

# logistic_aic <- extractAIC(logistic.rms)[2]
# logistic_aic # AIC = 6450.456


aic_mat <- matrix(c( log_log_aic, log_normal_aic, wei_aic, gomp_aic, expo_aic ), ncol = 1)

colnames(aic_mat) <- c('AIC')
rownames(aic_mat) <- c('Log-Logistic','Log-Normal','Weibull (AFT)','Gompertz','Exponential')
aic_mat



library(kableExtra)

# Convert matrix to a data frame
df_aic_mat <- as.data.frame(aic_mat)


caption.0 <- "<b>Akaike's information criterion for the parametric survival models fitted in the 50% Scenario</b>"

# Create a styled table using kableExtra
table.0 <- kable(df_aic_mat, format = "html", table.attr = "class='table'", 
               caption = caption.0) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Display the table
table.0











## Comparisons
# From the full data estimate the probability of survival up to 337 months (print the survival estimates per month)
observed <- summary(survfit(Surv(Time,event) ~ 1, data = data), times = seq(0,337,1)  )
head(observed)


# Extract the observed survival probabilities per month
observed_surv <- observed$surv
observed_surv <- observed_surv[161:337] # Select the observed survival probs after month 160
head(observed_surv)

summary(observed_surv)


# Select the survival estimates from month 161 up to 337
mean_surv_month_expo <- mean_surv_month_expo[1:177]


# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Exponential model fitted with RMS package
total.expo.rms <- as.matrix(observed_surv - mean_surv_month_expo) # Transformation needed for norm function to work


# Second way to calculate l1 norm via norm function
# l1.norm.expo.rms <- norm(total.expo.rms, type = "1") # specifies the one norm (maximum absolute column sum or l1 norm)
# l1.norm.expo.rms # 15.62566


# Frobenius_norm.expo.rms <- norm(total.expo.rms, type="F") # specifies the Frobenius norm (the Euclidean norm of total.expo treated as if it were a vector)
# Frobenius_norm.expo.rms # 1.49491


# Calculate l2_norm: The sum of the squared, observed minus the predicted survival probabilities for the Exponential model
l2_norm_expo.rms <- sqrt(sum(total.expo.rms^2))
l2_norm_expo.rms


# Calculate l1_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_expo.rms <- sum(abs(total.expo.rms))
l1_norm_expo.rms



# Calculate the mean Brier Score for the Exponential model fitted with RMS package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_expo.rms <- (1/n)*sum( (observed_surv - mean_surv_month_expo)^2 )
mean_squared_error_expo.rms


# Second way to calculate the mean BS for the Exponential model fitted with RMS package
mean_BS_expo.rms <- mean((observed_surv - mean_surv_month_expo)^2)
mean_BS_expo.rms









# Select the estimated survival probabilities from month 161 up to 337
expo_survivor_func <- expo_survivor_func[1:177]

# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Exponential model fitted with survHE package
total.expo.manual <- as.matrix(observed_surv - expo_survivor_func) # Transformation needed for norm function to work



# Calculate l2_norm: The sum of the squared, observed minus the predicted survival probabilities for the Exponential model
l2_norm_expo.manual <- sqrt(sum(total.expo.manual^2))
l2_norm_expo.manual


# Calculate l1_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_expo.manual <- sum(abs(total.expo.manual))
l1_norm_expo.manual



# Calculate the mean Brier Score for the Exponential model fitted with survHE package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_expo.manual <- (1/n)*sum( (observed_surv - expo_survivor_func)^2 )
mean_squared_error_expo.manual


# Second way to calculate the mean BS for the Exponential model fitted with survHE package
mean_BS_expo.manual <- mean((observed_surv - expo_survivor_func)^2)
mean_BS_expo.manual







# The fact that such a measure of overall accuracy of prediction may be decomposed into aspects of calibration and 
# discrimination has been particularly emphasized in a medical context
n <- length(expo_survivor_func)

# Unbiased = well calibrated
lack_of_calibration_expo_surv <- (1/n)*sum( (observed_surv - expo_survivor_func)*(1 - 2*expo_survivor_func) )
lack_of_calibration_expo_surv


# Sharpness = small variance
lack_of_spread_of_prediction_expo_surv <- (1/n)*sum( expo_survivor_func*(1 - expo_survivor_func ) )
lack_of_spread_of_prediction_expo_surv


# Expanded mean squared error
mean_squared_error_expanded_expo_surv <- lack_of_calibration_expo_surv^2 + lack_of_spread_of_prediction_expo_surv
mean_squared_error_expanded_expo_surv












# Extract the survival estimates from month 161 up to 337
mean_surv_month_wei <- mean_surv_month_wei.rms[1:177]


# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Weibull model fitted with RMS package
total.wei.rms <- as.matrix(observed_surv - mean_surv_month_wei) # Transformation needed for norm function to work

# One_norm.wei.rms <- norm(total.wei.rms, type = "1") # specifies the one norm (maximum absolute column sum)
# One_norm.wei.rms # same as l2_norm


# Calculate l1_norm: The sum of the squared, observed minus the predicted survival probabilities for the Exponential model
l2_norm_wei.rms <- sqrt(sum(total.wei.rms^2))
l2_norm_wei.rms 


# Calculate l2_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_wei.rms <- sum(abs(total.wei.rms))
l1_norm_wei.rms



# Calculate the mean Brier Score for the Weibull AFT model fitted with RMS package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_wei.rms <- (1/n)*sum( (observed_surv - mean_surv_month_wei)^2 )
mean_squared_error_wei.rms


# Second way to calculate the mean BS for the Weibull model fitted with RMS package
mean_BS_wei.rms <- mean((mean_surv_month_wei - observed_surv)^2)
mean_BS_wei.rms











# Extract the estimated survival probabilities from month 161 up to 337
weib_St <- weib$St[1:177]

# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Weibull model fitted with survHE package
total.wei.manual <- as.matrix(observed_surv - weib_St ) # Transformation needed for norm function to work


# Calculate l1_norm: The sum of the squared, observed minus the predicted survival probabilities for the Exponential model
l2_norm_wei.manual <- sqrt(sum(total.wei.manual^2))
l2_norm_wei.manual


# Calculate l2_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_wei.manual <- sum(abs(total.wei.manual))
l1_norm_wei.manual



# Calculate the mean Brier Score for the Weibull AFT model fitted with survival package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_wei.manual <- (1/n)*sum( (observed_surv - weib_St)^2 )
mean_squared_error_wei.manual


# Second way to calculate the mean BS for the Exponential model fitted with survHE package
mean_BS_wei.manual <- mean((weib_St - observed_surv)^2)
mean_BS_wei.manual








# Extract the estimated survival probabilities from month 161 up to 337
Log_Logistic_surv_St <- log.logis$St[1:177]

# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Log-Logistic model fitted with survival package
total.log.logis.manual <- as.matrix(observed_surv - Log_Logistic_surv_St ) # Transformation needed for norm function to work


# Calculate l1_norm: The sum of the squared, observed minus the predicted survival probabilities for the Log-Logistic model
l2_norm_log.logis.manual <- sqrt(sum(total.log.logis.manual^2))
l2_norm_log.logis.manual


# Calculate l2_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_log.logis.manual <- sum(abs(total.log.logis.manual))
l1_norm_log.logis.manual


# Calculate the mean Brier Score for the Log-Logistic model fitted with survival package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_log.logis.manual <- (1/n)*sum( (observed_surv - Log_Logistic_surv_St)^2 )
mean_squared_error_log.logis.manual












# Extract the estimated survival probabilities from month 161 up to 337
Log_Logistic_rms_St <- mean_surv_month_log.log.rms[1:177]

# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Log-Logistic model fitted with RMS package
total.log.logis.rms <- as.matrix(observed_surv - Log_Logistic_rms_St ) # Transformation needed for norm function to work


# Calculate l1_norm: The sum of the squared, observed minus the predicted survival probabilities for the Log-Logistic model
l2_norm_log.logis.rms <- sqrt(sum(total.log.logis.rms^2))
l2_norm_log.logis.rms


# Calculate l2_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_log.logis.rms <- sum(abs(total.log.logis.rms))
l1_norm_log.logis.rms



# Calculate the mean Brier Score for the Log-Logistic model fitted with RMS package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_log.logis.rms <- (1/n)*sum( (observed_surv - Log_Logistic_rms_St)^2 )
mean_squared_error_log.logis.rms















# Extract the estimated survival probabilities from month 161 up to 337
Log_Normal_rms_St <- mean_surv_month_log.norm.rms[1:177]

# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Log-Normal model fitted with RMS package
total.log.norm.rms <- as.matrix(observed_surv - Log_Normal_rms_St ) # Transformation needed for norm function to work


# Calculate l1_norm: The sum of the squared, observed minus the predicted survival probabilities for the Log-Logistic model
l2_norm_log.norm.rms <- sqrt(sum(total.log.norm.rms^2))
l2_norm_log.norm.rms


# Calculate l2_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_log.norm.rms <- sum(abs(total.log.norm.rms))
l1_norm_log.norm.rms



# Calculate the mean Brier Score for the Log-Normal model fitted with RMS package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_log.norm.rms <- (1/n)*sum( (observed_surv - Log_Normal_rms_St)^2 )
mean_squared_error_log.norm.rms







# Survivor function for the log normal model via the survHE package
x <- seq(161,400,1)
S.log.norm <- 1-pnorm(log(x), log_norm$coef[1] + log_norm$coef[2]*mean(age) + log_norm$coef[3]*mean(tumor_size) + log_norm$coef[4]*mean(lnep) +
                        log_norm$coef[5]*mean(gcs_2) + log_norm$coef[6]*mean(gcs_3) + log_norm$coef[7]*mean(gcs_4) +  
                        log_norm$coef[8]*mean(cancer_type_det_2) + log_norm$coef[9]*mean(cancer_type_det_3) + log_norm$coef[10]*mean(cancer_type_det_4) 
                      
                      , log_norm$scale)

# Extract the estimated survival probabilities from month 161 up to 337
Log_Normal_survHE_st <- S.log.norm[1:177] 


# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Log-Normal model fitted with RMS package
total.log.norm.survHE <- as.matrix(observed_surv - Log_Normal_survHE_st ) # Transformation needed for norm function to work


# Calculate l1_norm: The sum of the squared, observed minus the predicted survival probabilities for the Log-Logistic model
l2_norm_log.norm.survHE <- sqrt(sum(total.log.norm.survHE^2))
l2_norm_log.norm.survHE


# Calculate l2_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_log.norm.survHE <- sum(abs(total.log.norm.survHE))
l1_norm_log.norm.survHE



# Calculate the mean Brier Score for the Log-Normal model fitted with RMS package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_log.norm.survHE <- (1/n)*sum( (observed_surv - Log_Normal_survHE_st)^2 )
mean_squared_error_log.norm.survHE

















# Extract the estimated survival probabilities from month 161 up to 337
Gompertz_St <- gomp_survivor_func[1:177]

# One, Frobenius, l1 and l2 Norms per month regarding the Observed Survival Probs minus the Predicted Survival Probs obtained via the Gompertz model fitted with survHE package
total.gompertz.st <- as.matrix(observed_surv - Gompertz_St ) # Transformation needed for norm function to work


# Calculate l1_norm: The sum of the squared, observed minus the predicted survival probabilities for the Log-Logistic model
l2_norm_gmprtz_survHE <- sqrt(sum(total.gompertz.st^2))
l2_norm_gmprtz_survHE


# Calculate l2_norm: The sum of the absolute, observed minus the predicted survival probabilities for the Exponential model
l1_norm_gmprtz_survHE <- sum(abs(total.gompertz.st))
l1_norm_gmprtz_survHE



# Calculate the mean Brier Score for the Log-Normal model fitted with RMS package
# The mean squared error, also known in the context of probabilistic predictions as mean Brier Score or the mean probability score:
n <- length(observed_surv)

mean_squared_error_gmprtz_survHE <- (1/n)*sum( (observed_surv - Gompertz_St)^2 )
mean_squared_error_gmprtz_survHE














# Create a matrix to store the results of the various norms + BS score for the various models fitted
Norm_BS_mat <- matrix(round(c(l1_norm_expo.rms,l1_norm_expo.manual,l1_norm_wei.rms,l1_norm_wei.manual,l1_norm_log.logis.rms,  l1_norm_log.logis.manual, l1_norm_log.norm.rms, l1_norm_log.norm.survHE  ,l1_norm_gmprtz_survHE,
                              
                              l2_norm_expo.rms, l2_norm_expo.manual, l2_norm_wei.rms, l2_norm_wei.manual, l2_norm_log.logis.rms, l2_norm_log.logis.manual, l2_norm_log.norm.rms, l2_norm_log.norm.survHE  ,l2_norm_gmprtz_survHE,
                              mean_squared_error_expo.rms, mean_squared_error_expo.manual, mean_squared_error_wei.rms, mean_squared_error_wei.manual, 
                              mean_squared_error_log.logis.rms, mean_squared_error_log.logis.manual, mean_squared_error_log.norm.rms,  mean_squared_error_log.norm.survHE ,mean_squared_error_gmprtz_survHE
                              
),3),
ncol = 3, nrow = 9, byrow = F)

rownames(Norm_BS_mat) <- c('Expo (RMS)','Expo (survHE)','Weibull AFT (RMS)', 'Weibull AFT (survHE)', 'Log-Logistic (RMS)',
                           'Log-Logistic (survHE)','Log-Normal (RMS)','Log-Normal (survHE)', 'Gompertz (survHE)')
colnames(Norm_BS_mat) <- c('l1 Norm','l2 Norm','Mean Brier Score' )
Norm_BS_mat





library(kableExtra)

# Convert matrix to a data frame
df_Norm_BS_mat <- as.data.frame(Norm_BS_mat)


caption.1 <- "<b>Metrics to assess the predicted survival probabilities obtained by the parametric models for the 50% Scenario</b>"

# Create a styled table using kableExtra
table.1 <- kable(df_Norm_BS_mat, format = "html", table.attr = "class='table'", 
               caption = caption.1) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Display the table
table.1





















# Create a matrix to store the results of the various norms + BS score for the various models fitted
Norm_BS_mat_new <- matrix(round(c(l1_norm_expo.manual, l1_norm_wei.manual, l1_norm_log.logis.manual, l1_norm_log.norm.survHE  ,l1_norm_gmprtz_survHE,
                              
                              l2_norm_expo.manual, l2_norm_wei.manual, l2_norm_log.logis.manual,  l2_norm_log.norm.survHE  ,l2_norm_gmprtz_survHE,
                              mean_squared_error_expo.manual, mean_squared_error_wei.manual, 
                             mean_squared_error_log.logis.manual,  mean_squared_error_log.norm.survHE ,mean_squared_error_gmprtz_survHE
                              
),3),
ncol = 3, nrow = 5, byrow = F)

rownames(Norm_BS_mat_new) <- c('Exponential', 'Weibull AFT', 'Log-Logistic',
                           'Log-Normal', 'Gompertz')
colnames(Norm_BS_mat_new) <- c('l1 Norm','l2 Norm','Mean Brier Score' )
Norm_BS_mat_new





library(kableExtra)

# Convert matrix to a data frame
df_Norm_BS_mat_new <- as.data.frame(Norm_BS_mat_new)


caption.1.1 <- "<b>Metrics to assess the prediction for the 50% scenario</b>"

# Create a styled table using kableExtra
table.1.1 <- kable(df_Norm_BS_mat_new, format = "html", table.attr = "class='table'", 
                 caption = caption.1.1) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Display the table
table.1.1



##############################################      Mean Survival Time for the Exponential Model


# Closed form of the E(t) = MST for the Exponential model
1/lambda1



## summary
# t: the vector of times to be used in the computation. Default = NULL, which means
# the observed times will be used. NB: the vector of times should be: i) long
# enough so that S(t) goes to 0; and ii) dense enough so that the approximation to the AUC is sufficiently precise
# The summary method computes the mean survival time for the various models fitted -> To do: For the whole time horizon!)
MST_expo_fifty_percent_scenario <- summary(m1, mod = 1, t = seq(0,3300,1) )
MST_expo_fifty_percent_scenario <- MST_expo_fifty_percent_scenario$tab
MST_expo_fifty_percent_scenario




# Approximate the AUC via the Simpson's rule
library(Bolstad2)

lambda1
x <- seq(0,3300,1)
fx1 <- exp(-lambda1*x)

# Check that the survival probability at the last 50 months is near zero
tail(fx1,50)

estimate1 = sintegral(x,fx1)$int
estimate1



names(lambda1) <- NULL
data.frame( Estimated = 1/lambda1, Trapezoidal_rule_estimate = MST_expo_fifty_percent_scenario[1], Simpson_estimate = estimate1 )

mean_exp(lambda1)


##############################################      Mean Survival Time for the Weibull (AFT) Model



## summary
MST_wei_fifty_percent_scenario <- summary(m1, mod = 2, t = seq(0,3300,1))
MST_wei_fifty_percent_scenario <- MST_wei_fifty_percent_scenario$tab
MST_wei_fifty_percent_scenario


# Approximate the AUC via the Simpson's rule
library(Bolstad2)

x <- seq(0,3300,1)
fx2 <- exp(-(lambda33*x)^(1/p11))

  

# Check that the survival probability at the last 50 months is near zero
tail(fx2,50)

estimate2 = sintegral(x,fx2)$int
estimate2


data.frame( Trapezoidal_rule_estimate = MST_wei_fifty_percent_scenario[1], Simpson_estimate = estimate2 )







##############################################      Mean Survival Time for the Log-Logistic Model



## summary
MST_Log_Logistic_fifty_percent_scenario <- summary(m1, mod = 3, t = seq(0,3300,1))
MST_Log_Logistic_fifty_percent_scenario <- MST_Log_Logistic_fifty_percent_scenario$tab
MST_Log_Logistic_fifty_percent_scenario


# Approximate the AUC via the Simpson's rule
library(Bolstad2)


x <- seq(0,3300,1)
fx3 <- 1/( 1 + (lambda5*x)^scale )



# Check that the survival probability at the last 50 months is near zero
tail(fx3,50)

estimate3 = sintegral(x,fx3)$int
estimate3


data.frame( Trapezoidal_rule_estimate = MST_Log_Logistic_fifty_percent_scenario[1], Simpson_estimate = estimate3 )




##############################################      Mean Survival Time for the Log-Normal Model



## summary
MST_Log_Normal_fifty_percent_scenario <- summary( m1, mod = 4, t = seq(0,3300,1) )
MST_Log_Normal_fifty_percent_scenario <- MST_Log_Normal_fifty_percent_scenario$tab
MST_Log_Normal_fifty_percent_scenario







x <- seq(0,3300,1)
fx4 <- 1-pnorm(log(x), log_norm$coef[1] + log_norm$coef[2]*mean(age) + log_norm$coef[3]*mean(tumor_size) + log_norm$coef[4]*mean(lnep) +
log_norm$coef[5]*mean(gcs_2) + log_norm$coef[6]*mean(gcs_3) + log_norm$coef[7]*mean(gcs_4) +  
  log_norm$coef[8]*mean(cancer_type_det_2) + log_norm$coef[9]*mean(cancer_type_det_3) + log_norm$coef[10]*mean(cancer_type_det_4) 
  
  , log_norm$scale   )



# Check that the survival probability at the last 50 months is near zero
head(fx4,10)
tail(fx4,50)

estimate4 = sintegral(x,fx4)$int
estimate4











##############################################      Mean Survival Time for the Gompertz Model



## summary
MST_Gompertz_fifty_percent_scenario <- summary(m1, mod = 5, t = seq(0,3300,1))
MST_Gompertz_fifty_percent_scenario <- MST_Gompertz_fifty_percent_scenario$tab
MST_Gompertz_fifty_percent_scenario


# Approximate the AUC via the Simpson's rule
library(Bolstad2)


x <- seq(0,3300,1)
fx5 <- exp( -(b2/a2)*(exp(a2*x)-1))



# Check that the survival probability at the last 50 months is near zero
head(fx5,20)
tail(fx5,50)

estimate5 = sintegral(x,fx5)$int
estimate5


data.frame( Trapezoidal_rule_estimate = MST_Gompertz_fifty_percent_scenario[1], Simpson_estimate = estimate5 )






MST_fifty_percent_scenario_S <- rbind( Mean_Survival_Time_Exponential = estimate1,
                                     Mean_Survival_Time_Weibull = estimate2,
                                     Mean_Survival_Time_Log_Logistic = estimate3,
                                     Mean_Survival_Time_Log_Normal = estimate4,
                                     Mean_Survival_Time_Gompertz = estimate5
)

MST_fifty_percent_scenario_S



MST_fifty_percent_scenario_TR <- rbind( Mean_Survival_Time_Exponential = MST_expo_fifty_percent_scenario[1],
                                     Mean_Survival_Time_Weibull = MST_wei_fifty_percent_scenario[1],
                                     Mean_Survival_Time_Log_Logistic = MST_Log_Logistic_fifty_percent_scenario[1],
                                     Mean_Survival_Time_Log_Normal = MST_Log_Normal_fifty_percent_scenario[1],
                                     Mean_Survival_Time_Gompertz = MST_Gompertz_fifty_percent_scenario[1]
)

MST_fifty_percent_scenario_TR













is.matrix(MST_fifty_percent_scenario)


rownames(MST_fifty_percent_scenario) <- c("Exponential","Weibull (AFT)","Log-Logistic","Log-Normal","Gompertz")
colnames(MST_fifty_percent_scenario) <- "Mean Survival Time Estimate"
MST_fifty_percent_scenario




library(kableExtra)

# Convert matrix to a data frame
df_MST_fifty_percent_scenario <- as.data.frame(MST_fifty_percent_scenario)
df_MST_fifty_percent_scenario



# Create a fancy matrix with the mean survival time estimates for the 50% scenario
caption.2 <- "<b>Mean survival time estimates for the 50% scenario</b>"
# Create a styled table using kableExtra
table.2 <- kable(df_MST_fifty_percent_scenario, format = "html", table.attr = "class='table'", 
               caption = caption.2) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Display the table
table.2








########################################################          Restricted Mean Survival Time Estimates


# KM for the cut off data
km_fifty_percent_scenario <- survfit(Surv(Time, event) ~ 1, type = "kaplan-meier", data = cut_off_dat)

# Plot the KM while expanding the xlim to 400 months
plot(km_fifty_percent_scenario, conf.int = F, lwd = 2, xlab = "Months", main = 'Kaplan Meyer Plot', ylab = 'S(t)', xlim = c(0,400), las = 1 ) 


# Non-parametric estimate of the restricted mean survival time
print(km_fifty_percent_scenario, print.rmean = T)
rmean_fifty_percent_scenario <- 117




# Parametric Estimates of the restricted mean survival time
library(cubature)


# Define the function 
integrand_expo <- function(x) {
  S <- exp(-(lambda1 * x)) # The Exponential Survivor Function
  return(S)
}

# Evaluate the integral of the Exponential survivor function up to time 160
expo_estimate_rmean_fifty_percent_scenario <- adaptIntegrate(integrand_expo, lowerLimit = 0, upperLimit = 160)$integral
expo_estimate_rmean_fifty_percent_scenario

rmst_exp(t = 160,rate = lambda1) # another way, same result as above




# Define the function 
integrand_wei <- function(x) {
  S <- exp(-(lambda3*x)^pp) # The Weibull AFT Survivor Function
  return(S)
}

# Evaluate the integral of the Weibull AFT survivor function up to time 160
wei_estimate_rmean_fifty_percent_scenario <- adaptIntegrate(integrand_wei, lowerLimit = 0, upperLimit = 160)$integral
wei_estimate_rmean_fifty_percent_scenario







# Define the function 
integrand_log_logis <- function(x) {
  S <- 1/( 1 + (lambda5*x)^scale ) # The Log-Logistic Survivor Function
  return(S)
}

# Evaluate the integral of the Log-Logistic survivor function up to time 160
Log_logis_estimate_rmean_fifty_percent_scenario <- adaptIntegrate(integrand_log_logis, lowerLimit = 0, upperLimit = 160)$integral
Log_logis_estimate_rmean_fifty_percent_scenario







# Define the function 
integrand_gmprtz <- function(x) {
  S <-   exp( -(b2/a2)*(exp(a2*x)-1)) # The Gompertz Survivor Function
    return(S)
}

# Evaluate the integral of the Gompertz survivor function up to time 160
Gompertz_estimate_rmean_fifty_percent_scenario <- adaptIntegrate(integrand_gmprtz, lowerLimit = 0, upperLimit = 160)$integral
Gompertz_estimate_rmean_fifty_percent_scenario








# First way to calculate RMST for the Log-Normal model
print(log_norm)

s <- log_norm$scale
s

mu <- log_norm$coefficients[1] + log_norm$coefficients[2]*mean(age) + log_norm$coefficients[3]*mean(tumor_size) +
  log_norm$coefficients[4]*mean(lnep) + log_norm$coefficients[5]*mean(gcs_2) + log_norm$coefficients[6]*mean(gcs_3) +
  log_norm$coefficients[7]*mean(gcs_4) + log_norm$coefficients[8]*mean(cancer_type_det_2) + 
  log_norm$coefficients[9]*mean(cancer_type_det_3) + log_norm$coefficients[10]*mean(cancer_type_det_4)
mu


Log_Normal_rmean_fifty_percent_scenario <- rmst_lnorm(t = 160, mean = 5.096116, sdlog = 1.125707)
Log_Normal_rmean_fifty_percent_scenario





# Second way to calculate the RMST for the Log-Normal model
x <- seq(0,3300,1)
fx4 <- 1-pnorm(log(x), log_norm$coef[1] + log_norm$coef[2]*mean(age) + log_norm$coef[3]*mean(tumor_size) + log_norm$coef[4]*mean(lnep) +
                 log_norm$coef[5]*mean(gcs_2) + log_norm$coef[6]*mean(gcs_3) + log_norm$coef[7]*mean(gcs_4) +  
                 log_norm$coef[8]*mean(cancer_type_det_2) + log_norm$coef[9]*mean(cancer_type_det_3) + log_norm$coef[10]*mean(cancer_type_det_4) 
               
               , log_norm$scale) 

# Define the function 
integrand_log_norm <- function(x) {
  S <- 1-pnorm(log(x), log_norm$coef[1] + log_norm$coef[2]*mean(age) + log_norm$coef[3]*mean(tumor_size) + log_norm$coef[4]*mean(lnep) +
                 log_norm$coef[5]*mean(gcs_2) + log_norm$coef[6]*mean(gcs_3) + log_norm$coef[7]*mean(gcs_4) +  
                 log_norm$coef[8]*mean(cancer_type_det_2) + log_norm$coef[9]*mean(cancer_type_det_3) + log_norm$coef[10]*mean(cancer_type_det_4) 
               
               , log_norm$scale)  # The Log-Normal Survivor Function
  return(S)
}

# Evaluate the integral of the Gompertz survivor function up to time 160
Log_norm_estimate_rmean_fifty_percent_scenario <- adaptIntegrate(integrand_log_norm, lowerLimit = 0, upperLimit = 160)$integral
Log_norm_estimate_rmean_fifty_percent_scenario











# Store the restricted mean survival time estimates into a matrix
rmean_time_fifty_percent_scenario <- rbind( rmean_non_parametric_estimate = rmean_fifty_percent_scenario ,
                                            rmean_Exponential = expo_estimate_rmean_fifty_percent_scenario,
                                            rmean_Log_Normal = Log_Normal_rmean_fifty_percent_scenario,
                                            rMean_Log_Logistic = Log_logis_estimate_rmean_fifty_percent_scenario,
                                            rMean_Gompertz = Gompertz_estimate_rmean_fifty_percent_scenario,
                                                  rMean_Weibull = wei_estimate_rmean_fifty_percent_scenario
                                                 
                                                 
)






is.matrix(rmean_time_fifty_percent_scenario)


rownames(rmean_time_fifty_percent_scenario) <- c( "KM", "Exponential", "Log-Normal", "Log-Logistic", "Gompertz", "Weibull (AFT)" )
colnames(rmean_time_fifty_percent_scenario) <- "Restricted Mean Survival Time Estimate"
rmean_time_fifty_percent_scenario




library(kableExtra)

# Convert matrix to a data frame
df_rmean_time_fifty_percent_scenario <- as.data.frame(rmean_time_fifty_percent_scenario)
df_rmean_time_fifty_percent_scenario




library(knitr)
library(kableExtra)

caption.3 <- "<b>Restricted mean survival time estimates for the 50% scenario</b>"
table.3 <- kable(df_rmean_time_fifty_percent_scenario, format = "html", table.attr = "class='table'", 
                 caption = caption.3) %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE)

# Display the table
table.3











































