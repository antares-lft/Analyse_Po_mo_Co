# Chargement des packages
library(dplyr)
library(stringr)
library(tools)

# === 0. Définir les chemins ===
folder_path <- "/home/alafitte/Internship"

# Liste des espèces
species_list <- c(
  "Phyllostomus_discolor",
  "Phyllostomus_hastatus",
  "Myotis_davidii",
  "Myotis_myotis"
)

# Initialiser un dataframe vide
final_df <- data.frame(
  species = character(),
  pN_pS = numeric(),
  Pin = numeric(),
  Pis = numeric(),
  Nw_sum = numeric(),
  Sw_sum = numeric(),
  dN_dS = numeric(),
  dN = numeric(),
  dS = numeric(),
  Nw_sum.1 = numeric(),
  Sw_sum.1 = numeric(),
  taille_pop = numeric(),
  stringsAsFactors = FALSE
)

# === Boucle principale sur les espèces ===
for (sp in species_list) {
  
  # === 1. Lire pN/pS ===
  pnps_path <- file.path(folder_path, paste0(sp, "_pNpS_sim_4.txt"))
  pnps_lines <- readLines(pnps_path)
  pnps_lines <- pnps_lines[!grepl("^#", pnps_lines)]  # Ignore les commentaires
  
  if (length(pnps_lines) == 0) next
  
  last_pnps_line <- tail(pnps_lines, 1)
  pnps_cols <- strsplit(last_pnps_line, "\\s+")[[1]]
  pnps_values <- as.numeric(pnps_cols[1:5])  # Colonnes 1 à 5
  
  # === 2. Lire dN/dS ===
  dnds_path <- file.path(folder_path, paste0(sp, "_dNdS_sim_4.txt"))
  dnds_lines <- readLines(dnds_path)
  dnds_lines <- dnds_lines[!grepl("^#", dnds_lines)]
  
  if (length(dnds_lines) == 0) next
  
  last_dnds_line <- tail(dnds_lines, 1)
  dnds_cols <- strsplit(last_dnds_line, "\\s+")[[1]]
  dnds_values <- as.numeric(dnds_cols[1:5])  # Colonnes 1 à 5
  
  # === 3. Lire taille de population ===
  pop_path <- file.path(folder_path, paste0(sp, "_pop_sim_4.txt"))
  pop_lines <- readLines(pop_path)
  numeric_lines <- pop_lines[!grepl("^#", pop_lines)]
  
  if (length(numeric_lines) < 3) next
  
  pop_line3 <- numeric_lines[3]
  pop_cols <- strsplit(pop_line3, "\\s+")[[1]]
  taille_pop <- as.numeric(pop_cols[2])
  
  # === 4. Ajouter la ligne au tableau final ===
  final_df <- rbind(final_df, data.frame(
    species = sp,
    pN_pS = pnps_values[1],
    Pin = pnps_values[2],
    Pis = pnps_values[3],
    Nw_sum = pnps_values[4],
    Sw_sum = pnps_values[5],
    dN_dS = dnds_values[1],
    dN = dnds_values[2],
    dS = dnds_values[3],
    Nw_sum.1 = dnds_values[4],
    Sw_sum.1 = dnds_values[5],
    taille_pop = taille_pop,
    stringsAsFactors = FALSE
  ))
  
  # Affichage pour vérification
  cat("\n--- Lecture pour :", sp, "---\n")
  cat("pN/pS :", paste(pnps_values, collapse = ", "), "\n")
  cat("dN/dS :", paste(dnds_values, collapse = ", "), "\n")
  cat("Taille population :", taille_pop, "\n")
}

# === Résultat final ===
cat("\n=== Table finale ===\n")
print(final_df)

# Option : sauvegarde CSV
write.csv(final_df, file = "/home/alafitte/Internship/Rapport de stage/Données/resultats_final_4.csv", row.names = FALSE)
