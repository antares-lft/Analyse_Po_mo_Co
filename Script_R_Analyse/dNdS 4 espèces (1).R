library(dplyr)
library(ggplot2)

folder_path <- "/home/alafitte/Internship"
folder_path2 <- "/home/alafitte/Internship"

# Fichiers dNdS
file_dN_dS <- c(
  file.path(folder_path, "Phyllostomus_discolor_dNdS_sim_6.txt"),
  file.path(folder_path, "Phyllostomus_hastatus_dNdS_sim_6.txt"),
  file.path(folder_path, "Myotis_davidii_dNdS_sim_6.txt"),
  file.path(folder_path, "Myotis_myotis_dNdS_sim_6.txt")
)

# Fichiers pop
file_pop_sim <- c(
  file.path(folder_path2, "Phyllostomus_discolor_pop_sim_4.txt"),
  file.path(folder_path2, "Phyllostomus_hastatus_pop_sim_4.txt"),
  file.path(folder_path2, "Myotis_davidii_pop_sim_4.txt"),
  file.path(folder_path2, "Myotis_myotis_pop_sim_4.txt")
)

# Fonction pour lire et nettoyer dNdS
read_and_process <- function(file_path) {
  data <- read.table(file_path, header = TRUE, sep = "", stringsAsFactors = FALSE)
  data$dN.dS <- as.numeric(data$dN.dS)
  data <- data[is.finite(data$dN.dS), ]
  data$Species <- tools::file_path_sans_ext(basename(file_path))
  return(data)
}

# Fonction pour extraire la taille de population (ligne 4)
extract_number <- function(file_path) {
  lines <- readLines(file_path)
  line4 <- lines[4]
  parts <- unlist(strsplit(line4, "\\s+"))
  number <- as.numeric(parts[2])  # "p1 500 H" → on prend "500"
  return(number)
}

# Lecture des données dNdS
list_data <- lapply(file_dN_dS, read_and_process)
all_data <- bind_rows(list_data)

# Extraction des tailles de population
pop_sizes <- sapply(file_pop_sim, extract_number)

# Associer chaque espèce à sa pop
species_pop_df <- data.frame(
  Species = unique(all_data$Species),
  pop = pop_sizes
)

# Fusion avec les tailles de pop
all_data <- merge(all_data, species_pop_df, by = "Species")

# Moyennes et écart-type par population

summary_df <- all_data %>%
  group_by(pop) %>%
  summarise(
    mean_dnds = mean(dN.dS),
    sd_dnds = sd(dN.dS),
    n = n(),
    se = sd_dnds / sqrt(n),
    t_value = qt(0.975, df = n - 1),
    ci_lower = mean_dnds - t_value * se,
    ci_upper = mean_dnds + t_value * se)

# Tracer
ggplot(all_data, aes(x = factor(pop), y = dN.dS)) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "red") +
  geom_point(data = summary_df, aes(x = factor(pop), y = mean_dnds),
             color = "blue", size = 4, shape = 18) +
  geom_errorbar(data = summary_df,
                aes(x = factor(pop), ymin = ci_lower, ymax = ci_upper),
                width = 0.2, color = "black",
                inherit.aes = FALSE) +
  theme_minimal() +
  xlab("Population size") +
  ylab("dN/dS") +
  ggtitle("dN/dS distribution by population size with r=0.5 and mu=1e-6") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
