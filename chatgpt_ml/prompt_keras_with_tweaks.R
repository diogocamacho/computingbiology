# Load required libraries
library(tidyverse)
library(tidymodels)
library(keras3)

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
iris_split <- initial_split(iris, prop = 3/4)
iris_train <- training(iris_split)
iris_test <- testing(iris_split)

# Convert Species to numeric
iris_train$Species <- as.numeric(iris_train$Species) - 1
iris_test$Species <- as.numeric(iris_test$Species) - 1

# Define model
model <- keras_model_sequential() %>%
  layer_dense(units = 112/2, activation = 'relu', input_shape = ncol(iris_train[, 1:4]), name = "input") %>%
  layer_dense(units = 112/4, activation = 'relu', name = "deep1") %>%
  layer_dense(units = 112/8, activation = 'relu', name = "deep2") %>%
  # layer_dense(units = 4, activation = 'relu', name = "deep3") %>%
  layer_dense(units = 3, activation = 'softmax', name = "output")

# Compile model
model %>% compile(
  loss = 'sparse_categorical_crossentropy',
  optimizer = optimizer_adam(),
  metrics = c('accuracy')
)

# Train model
history <- model %>% 
  fit(
  as.matrix(iris_train[,1:4]), iris_train$Species,
  epochs = 1000, 
  batch_size = nrow(iris_train) * 0.3,
  validation_split = 0.2,
  callbacks = callback_early_stopping(monitor = "val_loss", 
                                      restore_best_weights = TRUE, 
                                      start_from_epoch = 1,  
                                      patience = 10, verbose = 1, mode = "auto")
  )

# Evaluate model
model %>% evaluate(as.matrix(iris_test[,1:4]), iris_test$Species)
