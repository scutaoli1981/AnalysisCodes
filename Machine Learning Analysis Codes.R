# Load necessary libraries
library(randomForest)
library(caret)
# Assuming 'data' is your dataset and 'target' is the response variable.
# Replace with actual variable names.
data <- read.csv('your_data.csv') # Replace with the path to your data file
target <- 'target_variable_name'  # Replace with your target variable name
# Prepare data for training and testing
data$target <- as.factor(data$target) # Convert target variable to factor if it's categorical
# Define control using 10-fold CV
fitControl <- trainControl(
  method = 'cv',
  number = 10,
  summaryFunction = twoClassSummary,
  classProbs = TRUE, # Set to TRUE if it's a classification problem
  verboseIter = TRUE
)
# Train the model using cross-validation
set.seed(123) # For reproducibility
rf_model <- train(
  formula(paste(target, '~ .')), # Target variable and predictors
  data = data,
  method = 'rf', # Method for Random Forest
  trControl = fitControl,
  metric = 'Accuracy' # Change this if you have a different performance metric in mind
)
# Output the results
print(rf_model)
#Load necessary library
library(pROC)
#Assuming you have an actual response vector and a prediction score or probability vector response that are real classification outcomes, predictions are probability examples of model predictions.
response <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,1,0,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0)
predictions <- c(0.7,0.4,0.9,0.3,0.1,0.4,0.3,0.3,0.2,0.4,0.3,0.1,0.2,0.4,0,0.3,0.2,0.2,0.8,0.7,0,0.6,0,0.7,0.4,0.4,0,0.2,0.4,0,0.7,0.4,0.7,0.2,0.8,0.1,0.6,0.8,0.3,0.4,0.7,0.9,0.2,0.8,0.2,0.4,0.2,0.1,0.8,0.7,0.9,0.4,0.8,0.1,0,0.4,0.1,0.2,0.1,0.2)
#Create ROC curve objects
roc_obj <- roc(response, predictions)
#Calculate AUC
auc_value <- auc(roc_obj)
#Print AUC value
print(auc_value)
#Plot ROC curve
plot(roc_obj, main="ROC Curve", col="#1c61b6")
