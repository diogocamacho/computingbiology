# Load required libraries
library(tidyverse)
library(tidymodels)
library(caret)

# Load iris dataset
data(iris)

# EDA
summary(iris)
pairs(iris)

# PCA
iris.pca <- prcomp(iris[,1:4], center = TRUE, scale. = TRUE)
summary(iris.pca)
plot(iris.pca, type = "l")

# Data preprocessing
iris_split <- initial_split(iris, prop = 0.75)
iris_train <- training(iris_split)
iris_test <- testing(iris_split)

# Define cross-validation method
cv <- trainControl(method = "cv", number = 10)

# Define models
models <- c("rf", "svmRadial", "knn", "nnet", "rpart")

# Initialize a list to store results
results <- list()

# Train models and store results
for(model in models) {
  set.seed(123)
  fit <- train(Species ~ ., data = iris_train, method = model, trControl = cv)
  pred <- predict(fit, newdata = iris_test)
  cm <- confusionMatrix(pred, iris_test$Species)
  results[[model]] <- list(fit = fit, cm = cm)
}

# Print results
for(model in models) {
  print(paste("Model:", model))
  print(results[[model]]$cm)
}
