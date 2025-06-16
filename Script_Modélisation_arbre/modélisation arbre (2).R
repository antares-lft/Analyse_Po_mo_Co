library(ape)
library(ggplot2)
library(ggtree)

# Lire le fichier Newick
tree <- read.tree("votre_fichier.newick")

# Exemple de fonction pour extraire les effectifs de population des noms des nœuds
# Supposons que le nom du nœud est sous la forme "Nom_Effectif"
extract_population_size <- function(node_name) {
  parts <- strsplit(node_name, "_")[[1]]
  if (length(parts) > 1) {
    as.numeric(parts[2])
  } else {
    NA
  }
}

# Appliquer la fonction à tous les nœuds de l'arbre
population_sizes <- sapply(tree$tip.label, extract_population_size)

# Créer un data frame avec les informations des branches
branches <- data.frame(
  branch = 1:length(tree$edge),
  edge = tree$edge,
  length = tree$edge.length,
  population_size = NA
)

# Associer les effectifs de population aux branches
# Cela dépend de la structure de votre arbre et de la manière dont les effectifs sont associés aux branches
# Vous devrez adapter cette partie en fonction de votre fichier Newick
branches$population_size <- population_sizes[branches$edge[, 2]]

# Créer un objet ggtree
p <- ggtree(tree, size=0.5) + theme_tree2()

# Ajouter les effectifs de population comme épaisseur des branches
p <- p + aes(size=population_size) +
  scale_size(range = c(0.5, 5)) + # Ajustez la plage de taille selon vos besoins
  geom_tree()

# Afficher l'arbre
print(p)