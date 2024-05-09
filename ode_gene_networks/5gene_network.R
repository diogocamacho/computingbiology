loadSBML("~/Desktop/tmp_gene_network_copasi.xml")

pars <- getParameters()$key
pars <- pars[grep("V", pars)]
sp <- setdiff(getSpecies()$name, "species")

basal_state <- CoRC::runSteadyState()$species

##### PERTURBATIONS: KO OF EACH GENE #####
v_value <- 0

D <- matrix(0, nrow = length(sp) * length(v_value), ncol = length(sp))

for(i in seq(1, nrow(D))) {
  for (k in seq(1, length(v_value))) {
    CoRC::setParameters(key = pars[i], value = v_value[k])
    kod <- runSteadyState()$species
    D[i, which(sp %in% kod$name)] <- kod$concentration
    
  }
  CoRC::setParameters(key = pars[i], value = 1)
}
D[D < 1e-5] <- 0
rownames(D) <- c("G1_ko", "G2_ko", "G3_ko", "G4_ko", "G5_ko")
colnames(D) <- c("G1", "G2", "G3", "G4", "G5")


##### PERTURBATIONS: KO OF PAIRS #####
all_pairs <- expand.grid(pars, pars)
all_pairs <- all_pairs[all_pairs$Var1 != all_pairs$Var2, ]
v_value <- 0


Dpairs <- matrix(0, nrow = nrow(all_pairs) * length(v_value), ncol = length(sp))

for(i in seq(1, nrow(Dpairs))) {
  for (k in seq(1, length(v_value))) {
    CoRC::setParameters(key = as.character(all_pairs[i, 1]), value = v_value[k])
    CoRC::setParameters(key = as.character(all_pairs[i, 2]), value = v_value[k])
    kod <- runSteadyState()$species
    Dpairs[i, which(sp %in% kod$name)] <- kod$concentration
  }
  CoRC::setParameters(key = as.character(all_pairs[i, 1]), value = 1)
  CoRC::setParameters(key = as.character(all_pairs[i, 2]), value = 1)
  
}



##### PERTURBATIONS: ONE NODE, RANDOM #####
drug_concs <- seq(0, .5, by = 0.01)
pert_exp <- vector(mode = "list", length = 1000)
experimental_noise <- seq(0.001, 0.05, by = 0.001)
experimental_noise <- c((1+ experimental_noise), (1 - experimental_noise))

for(i in seq(1, 1000)) {
  
  dc <- sample(drug_concs, 1)
  drug_target <- sample(pars, 1)
  noisy_v <- sample(experimental_noise, 5)
  CoRC::setParameters(key = pars, value = noisy_v)
  CoRC::setParameters(key = drug_target, value = dc)
  
  dexp <- runSteadyState()$species
  
  pert_exp[[i]] <- tibble::tibble(perturbed_parameter = drug_target,
                               perturbed_species = substring(text = drug_target, first = 2, last = 3),
                               parameter_value = dc,
                               species = dexp$name,
                               steady_state = dexp$concentration)

}
pert_exp <- do.call(rbind, pert_exp)

# diagnostic plots
pert_exp %>% 
  dplyr::count(., perturbed_species) %>% 
  ggplot(aes(x = perturbed_species, y = n, fill = perturbed_species)) +
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_viridis_d() + 
  labs(x = "Gene perturbed", y = "Number of times perturbed") + 
  theme_bw() + 
  theme(panel.background = element_blank())

# example gene: G2
pert_exp %>% 
  dplyr::filter(species == "G2", perturbed_species != "G2") %>%
  ggplot(aes(x = perturbed_species, y = steady_state, fill = perturbed_species)) + 
  geom_boxplot() + 
  labs(x = "Perturbation", y = "G2 steady state") + 
  theme_bw() + 
  theme(panel.background = element_blank())

pert_exp %>%
  dplyr::filter(., perturbed_species == "G1", species != "G1") %>%
  ggplot(aes(x = parameter_value, y = steady_state, fill = species)) + 
  geom_point(shape = 21, color = "black", size = 2, alpha = 0.75) + 
  facet_wrap(. ~ species, scales = "free") + 
  labs(x = "V_G1 parameter", y = "Species steady state") + 
  theme_bw() +
  theme(panel.background = element_blank())


# reset model
CoRC::setParameters(key = pars, value = 1)
