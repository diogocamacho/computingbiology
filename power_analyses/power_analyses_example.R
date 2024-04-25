library(RnaSeqSampleSize)


# how many samples do i need?
RnaSeqSampleSize::sample_size(power = 0.8, m = 12148, m1 = 197, f = 0.05, k = 1, )  

# let's map samples vs power
power_level <- seq(0.01, 0.99, by = 0.01)
num_samples <- vector(length = length(power_level))
for (i in seq(1, length(power_level))) {
  num_samples[i] <- RnaSeqSampleSize::sample_size(power = power_level[i], m = 12148, m1 = 197, f = 0.05, k = 1, )  
}

# build data frame
power_analyses <- tibble::tibble(power_level = power_level,
                                 number_samples = num_samples)

power_analyses %>% 
  ggplot(aes(x = number_samples, y = power_level*100)) +
  geom_point(shape = 21, fill = "red", size = 2, alpha = 0.75) +
  labs(x = "Number of Samples", y = "Power (%)") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black"))