# Charge packages
library(treeio)
library(ggtree)

# Read the tree
arbre <- read.nhx("/home/alafitte/Internship/Rapport de stage/Arbre/bat_genes_complete_filtered-PhyML_tree.nhx")

# Show data associed to each node
head(arbre@data)

# Visualize the tree
p <- ggtree(arbre, aes(color = as.numeric(W))) + 
  geom_tiplab() + 
  scale_color_viridis_c(option = "plasma") + 
  theme_tree2() + 
  ggtitle("Visualisation phylogenetic tree")

# Showing the tree
print(p)
