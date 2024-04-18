set.seed(08030727)

# set our "truth" up
# we're assuming that we can measure 20000 things, and that our true signal comes from ~10% of them
true_rate <- 0.1
#rbinom(n = 20000, size = 1, prob = true_rate)
rnbinom(n = 20000, size = 1, prob = true_rate)
samp_vars <- sample(x = seq(1, 20000), size = length(seq(1,20000)) * true_rate, replace = FALSE)
truth_vector <- vector(length = length(seq(1,20000)))
truth_vector[samp_vars] <- 1

# now let's explore accuracy of a model as we change the number of things perturbed, under the same distribution
p <- seq(0.05, 0.95, by = 0.05)

# let's establish an experimental error rate: something that that you'd expect to see based on experimental error alone, independent of model
# for this example, we'll set that at 5%
err <- 0.05
  
# this is just a variable where we will store our results
ntimes <- 1000 # how many times we will do this
model_performance <- vector(mode = "list", length = ntimes)


# now let's do some simulations
for (i in seq(1, ntimes)) {
  
  # where can i sample from
  tv <- which(truth_vector == 1)
  
  # sample how much of the truth vector you will cover
  s <- sample(p, size = 1)
  
  # sample the truth vector based on s
  this_sampling <- sample(tv, size = length(which(truth_vector == 1)) * s, replace = FALSE)
  
  # add additional observations to complete length of truth vector observations
  add_obs <- sample(x = setdiff(seq(1, length(truth_vector)), this_sampling), 
                    size = length(tv) - length(this_sampling), 
                    replace = FALSE)
  
  # now we will sample, from a binomial, some random error --> our experimental error
  exp_err <- sample(x = seq(1, 20000), size = length(seq(1,20000)) * err, replace = FALSE)
  
  # now we put these things together
  obs <- vector(mode = "integer", length = 20000)
  obs[this_sampling] <- 1 # <-- these are the true things we measured
  obs[add_obs] <- 1 # <-- these are the other things we measured (our false positives)
  obs <- obs + obs[exp_err] # <-- this is our complete "perturbation"
  obs[obs > 1] <- 1 # a gimmick to make things simple
  
  
  tp <- length(intersect(which(x == 1), which(obs == 1)))
  tn <- length(intersect(which(x == 0), which(obs == 0)))
  fp <- length(intersect(which(x == 0), which(obs == 1)))
  fn <- length(intersect(which(x == 1), which(obs == 0)))
  
  model_performance[[i]] <- tibble::tibble(truth_percentage = s,
                                           number_hits = length(which(obs == 1)),
                                           tp = tp,
                                           tn = tn,
                                           fp = fp,
                                           fn = fn,
                                           acc = (tp + tn) / (tp + tn + fp + fn),
                                           ppv = tp / (tp + fp), # precision
                                           fdr = 1 - ppv,
                                           tpr = tp / (tp + fn), # recall, sensitivity
                                           tnr = tn / (tn + fp), # specificity
                                           fnr = fn / (tp + fn), # type II errors
                                           fpr = fp / (fp + tn), # type I error
                                           f1score = (2 * ppv * tpr) / (ppv + tpr)
                                           )
  
}
model_performance <- do.call(rbind, model_performance)

# let's plot some things
# does our model capture our error consistently?
model_performance %>% 
  ggplot(aes(x = truth_percentage, y = fpr)) + 
  geom_point(shape = 19)

# how does our sensitivity and specificity compare as we get closer to X?
model_performance %>% 
  ggplot(aes(x = tpr, y = ppv)) +
  geom_point(shape = 21)

# how does the F1 score change as we get closer to X?
model_performance %>% 
  ggplot(aes(x = truth_percentage, y = f1score)) +
  geom_point(shape = 21)

# how many things do we find compared to our truth?
model_performance %>% 
  ggplot(aes(x = as.factor(truth_percentage), y = tpr)) +
  geom_boxplot()



