# Show me how biology is hard


valproic50 <- which(cauldron::cmap_instances$name == "valproic acid" & 
      cauldron::cmap_instances$cell_line == "MCF7" & 
      cauldron::cmap_instances$concentration == 50 & 
      cauldron::cmap_instances$vehicle == "DMSO")


A <- cauldron::cmap_dataset$M[, valproic50]

number_pgenes <- apply(A, 2, function(x) length(which(x != 0)))
genes_pert <- apply(A, 1, function(x) length(which(x != 0)))