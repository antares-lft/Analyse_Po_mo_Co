# Installer les packages nécessaires
if (!requireNamespace("treeio", quietly = TRUE)) install.packages("treeio")
if (!requireNamespace("ggtree", quietly = TRUE)) install.packages("ggtree")

# Charger les packages
library(treeio)
library(ggtree)

# Lire l'arbre avec les annotations NHX
arbre <- read.nhx("/home/alafitte/Internship/Rapport de stage/Arbre/bat_genes_complete_filtered-PhyML_tree.nhx")

# Afficher la table des données associées à chaque nœud
head(arbre@data)

# Visualiser l'arbre avec ggtree et couleur selon la taille de population W
p <- ggtree(arbre, aes(color = as.numeric(W))) + 
  geom_tiplab() + 
  scale_color_viridis_c(option = "plasma") + 
  theme_tree2() + 
  ggtitle("Visualisation phylogenetic tree")

# Afficher
print(p)