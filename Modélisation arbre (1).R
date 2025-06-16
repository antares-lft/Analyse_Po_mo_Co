library(ggtree)
library(dplyr)
library(stringr)
library(scales)

# Lecture de l'arbre Newick (en texte ici, adapte selon fichier)
tree_file <- "/home/alafitte/Internship/Rapport de stage/Arbre/bat_genes_complete_filtered-PhyML_tree.nhx"
tree <- read.nhx(tree_file)

# 2. Extraire la taille de population W des labels (si pas déjà dans tree@data)
if (!"W" %in% colnames(tree@data)) {
  tree@data <- tree@data %>%
    mutate(
      W = as.numeric(str_extract(label, "(?<=W=)\\d+"))
    )
}

# 3. Remplacer NA par la moyenne pour éviter erreurs (optionnel : peut aussi choisir min ou autre)
tree@data$W[is.na(tree@data$W)] <- mean(tree@data$W, na.rm = TRUE)

# 4. Calculer une version log-scaled et mise à l’échelle de W pour la taille des branches
tree@data <- tree@data %>%
  mutate(
    W_scaled = rescale(log10(W + 1), to = c(0.5, 5))
  )

# 5. Tracer l'arbre avec épaisseur des branches selon W_scaled
p <- ggtree(tree, aes(size = W_scaled)) +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_tree2()

# 6. Ajouter les labels W sur les noeuds internes
p <- p + geom_text2(aes(subset = !isTip, label = W), hjust = -0.3)

# Afficher le graphique
print(p)