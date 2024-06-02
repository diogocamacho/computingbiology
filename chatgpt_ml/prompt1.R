# Load required libraries
library(tidyverse)
library(tidymodels)
library(caret)
library(xgboost)

# Load iris dataset
data(iris)

# Exploratory Data Analysis
# Summary of the dataset
summary(iris)

# Pair plot
pairs(iris[1:4], main = "Iris Data (red=setosa, green=versicolor, blue=virginica)",
      pch = 21, bg = c("red", "green3", "blue")[unclass(iris$Species)])

# Boxplot for each measurement
par(mfrow=c(2,2))
for(i in 1:4) {
  boxplot(iris[,i]~iris$Species, main=names(iris)[i])
}

# Data preprocessing
# Split the data into training and testing sets
set.seed(123)
iris_split <- initial_split(iris, prop = 0.75, strata = Species)
iris_train <- training(iris_split)
iris_test <- testing(iris_split)

# Recipe for preprocessing
iris_recipe <- recipe(Species ~ ., data = iris_train) %>%
  step_normalize(all_predictors())

# Model specification
xgb_spec <- boost_tree(trees = 1000, tree_depth = 5) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# Workflow
iris_workflow <- workflow() %>%
  add_recipe(iris_recipe) %>%
  add_model(xgb_spec)

# Model fitting
iris_fit <- fit(iris_workflow, data = iris_train)

# Model evaluation
iris_test_predictions <- predict(iris_fit, iris_test) %>%
  bind_cols(iris_test)
iris_metrics <- metrics(iris_test_predictions, truth = Species, estimate = .pred_class)

# Print out the results
print(iris_metrics)
