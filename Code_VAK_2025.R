# =============================================================================
# ECOFOR – Floristics analyses
# Manuscript: "Tree diversity of Atlantic Forest remnants: connectivity and size
#             as drivers of floristic richness" (Ramos et al.)
#
# This script reproduces core floristic community analyses:
# - community matrices (plot and subplot)
# - rarefaction/extrapolation (iNEXT)
# - Shannon diversity + Tukey contrasts (emmeans)
# - shared vs exclusive species (binary matrix + chord diagram)
# - NMDS (vegan) + convex hulls
# - UPGMA clustering (Bray-Curtis)
# - ANOSIM
# - phytosociological descriptors using fitoR (sourced from R/fitoR.R)
# =============================================================================

set.seed(1313)
rm(list = ls())

# ---- Paths ----
# Project structure suggested:
# .
# ├─ data/Dados.csv
# ├─ R/fitoR.R
# └─ analysis/01_ecofor_floristics.R
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

data_file <- here::here("data", "Dados.csv")

# ---- Packages ----
required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "readr", "stringr", "forcats",
  "vegan", "iNEXT", "emmeans", "multcomp",
  "viridis", "car", "caret",
  "ggdendro", "dendextend",
  "indicspecies",
  "chorddiag", "circlize",
  "knitr", "kableExtra", "scales"
)

missing_pkgs <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# ---- Load data ----
# NOTE: original code used read.csv2 with dec=",".
# readr::read_csv2 uses ';' separator by default and decimal ','.
dados <- readr::read_csv2(data_file, show_col_types = FALSE)

# ---- Expand PAP to PAP1, PAP2, ... per tree tag (Placa) ----
expand_pap_columns <- function(df, id_col = "Placa", value_col = "PAP") {
  df %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(PAP_index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(
      names_from = PAP_index,
      values_from = dplyr::all_of(value_col),
      names_prefix = "PAP"
    )
}

# original subset: dados[, c(4, 7)] (risky). Use names if possible:
# adjust these names if your columns differ.
data_expanded <- expand_pap_columns(dados %>% dplyr::select(Placa, PAP))

dados <- dados %>%
  dplyr::left_join(data_expanded, by = "Placa") %>%
  dplyr::distinct(Placa, .keep_all = TRUE)

# ---- Subplot id ----
dados <- dados %>%
  dplyr::mutate(subplot = paste0(Frag, "_", Sub))

# ---- Data cleaning ----
dados_clean <- dados %>%
  dplyr::filter(!familia %in% c("Morta", "Indeterminada"))

# ---- Abundance by spp x subplot (helper table) ----
abd_com <- dados_clean %>%
  dplyr::group_by(spp, subplot) %>%
  dplyr::summarise(
    abd_mean = mean(quantity, na.rm = TRUE),
    max_n = dplyr::n(),
    .groups = "drop"
  )

# ---- Community matrix: plot (Frag) x spp ----
comp_plot <- dados_clean %>%
  dplyr::group_by(Frag, spp) %>%
  dplyr::summarise(quantity = sum(quantity, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = spp, values_from = quantity, values_fill = 0) %>%
  tibble::column_to_rownames("Frag")

# ---- Remove singletons (species occurring once in the dataset) ----
sp_selection <- abd_com %>%
  dplyr::filter(max_n > 1) %>%
  dplyr::pull(spp) %>%
  unique()

dados_nosingle <- dados_clean %>%
  dplyr::filter(spp %in% sp_selection)

# ---- Community matrix: subplot x spp (full + no singletons) ----
comp_sub <- dados_clean %>%
  dplyr::group_by(subplot, spp) %>%
  dplyr::summarise(quantity = sum(quantity, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = spp, values_from = quantity, values_fill = 0) %>%
  tibble::column_to_rownames("subplot")

comp_clean_sub <- dados_nosingle %>%
  dplyr::group_by(subplot, spp) %>%
  dplyr::summarise(quantity = sum(quantity, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = spp, values_from = quantity, values_fill = 0) %>%
  tibble::column_to_rownames("subplot")

# =============================================================================
# Analyses
# =============================================================================

# ---- Rarefaction / extrapolation (q = 0) ----
rare <- iNEXT::iNEXT(t(comp_plot), q = 0, datatype = "abundance")
rare_result <- rare$iNextEst$size_based %>%
  dplyr::mutate(plot = rep(rownames(comp_plot), each = dplyr::n() / nrow(comp_plot)))

ggplot(rare_result, aes(x = m, y = qD, colour = plot, fill = plot)) +
  geom_line() +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.2, colour = NA) +
  theme_classic() +
  labs(x = "Sample size", y = "Species richness (q = 0)")

# ---- Shannon diversity (subplot level) ----
shannon_result <- data.frame(sha = vegan::diversity(comp_sub, index = "shannon")) %>%
  tibble::rownames_to_column("subplot") %>%
  dplyr::mutate(plot = stringr::str_remove(subplot, "_.*$"))  # assumes subplot = "Frag_Sub"

mod <- lm(sha ~ plot, data = shannon_result)

mod_means_contr <- emmeans::emmeans(mod, pairwise ~ plot, adjust = "tukey")
mod_letters <- multcomp::cld(mod_means_contr$emmeans, Letters = letters)
print(mod_letters)

ggplot(shannon_result, aes(x = plot, y = sha, colour = plot, fill = plot)) +
  geom_boxplot(alpha = 0.2) +
  theme_minimal() +
  labs(x = NULL, y = "Shannon diversity (H')") +
  guides(colour = "none", fill = "none")

# ---- Shared vs exclusive species (binary) + chord diagram ----
community_bin <- comp_plot
community_bin[community_bin > 0] <- 1

# Shared species between plots = crossproduct of binary presence/absence
shared_species_matrix <- as.matrix(community_bin) %*% t(as.matrix(community_bin))
diag(shared_species_matrix) <- 0

# Simple chord diagram (interactive htmlwidget)
chorddiag::chorddiag(shared_species_matrix)

# Exclusive species per plot = present only in that plot
exclusive_species <- rowSums(community_bin == 1 & colSums(community_bin) == 1)
total_species <- rowSums(community_bin > 0)
shared_species <- total_species - exclusive_species

species_data <- data.frame(
  Area = rownames(community_bin),
  Total = total_species,
  Exclusive = exclusive_species,
  Shared = shared_species
)

species_data_long <- species_data %>%
  tidyr::pivot_longer(cols = c("Exclusive", "Shared"), names_to = "Species_type", values_to = "Count")

ggplot(species_data_long, aes(x = Area, y = Count, fill = Species_type)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(Exclusive = "tomato4", Shared = "green4")) +
  labs(
    title = "Species richness by area: exclusive vs shared",
    x = "Area",
    y = "Number of species",
    fill = "Type"
  )

# ---- NMDS (no singletons; subplot level) ----
com_mds <- vegan::metaMDS(comp_clean_sub, k = 4, trymax = 100)
print(com_mds$stress)

sub_scores <- as.data.frame(vegan::scores(com_mds, display = "sites"))
sub_scores$subplot <- rownames(sub_scores)
sub_scores$plot <- stringr::str_remove(sub_scores$subplot, "_.*$")

create_hull <- function(df) {
  df[chull(df[, c("NMDS1", "NMDS2")]), , drop = FALSE]
}

hull_data <- sub_scores %>%
  dplyr::group_by(plot) %>%
  dplyr::group_modify(~create_hull(.x)) %>%
  dplyr::ungroup()

ggplot(sub_scores, aes(x = NMDS1, y = NMDS2, colour = plot)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_polygon(
    data = hull_data,
    aes(fill = plot, group = plot),
    alpha = 0.25,
    colour = NA
  ) +
  theme_minimal() +
  labs(title = "NMDS (Bray-Curtis)", x = "NMDS1", y = "NMDS2")

# ---- UPGMA clustering (Bray-Curtis) ----
dissimilarity_matrix <- vegan::vegdist(comp_plot, method = "bray")
upgma_clustering <- hclust(dissimilarity_matrix, method = "average")

dendro <- ggdendro::dendro_data(as.dendrogram(upgma_clustering))
segment_data <- ggdendro::segment(dendro)
label_data <- ggdendro::label(dendro)

ggplot() +
  geom_segment(data = segment_data,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label_data,
            aes(x = x, y = y, label = label),
            hjust = 1, size = 3) +
  coord_flip() +
  scale_y_reverse() +
  theme_minimal() +
  labs(title = "UPGMA dendrogram (Bray-Curtis)", x = "Areas", y = "Dissimilarity")

# ---- ANOSIM (subplot groups by plot) ----
group_vector <- shannon_result$plot[match(rownames(comp_clean_sub), shannon_result$subplot)]
anosim_res <- vegan::anosim(comp_clean_sub, group_vector)
print(anosim_res)
plot(anosim_res)

# =============================================================================
# Phytosociological descriptors (fitoR)
# =============================================================================

source(here::here("R", "fitoR.R"))

# Per plot (species-level)
dados_fito <- dados_clean %>%
  dplyr::select(Frag, Sub, spp, dplyr::starts_with("PAP")) %>%
  dplyr::mutate(parc = paste0(Frag, "_", Sub))  # required name in fitoR

# Convert circumference (PAP) to diameter (DAP): DAP = CAP / pi
# Here PAP seems to be CAP (circumference). Adjust if needed.
pap_cols <- grep("^PAP\\d+$", colnames(dados_fito), value = TRUE)
dap_mat <- as.data.frame(dados_fito[, pap_cols, drop = FALSE] / pi)

dados_fito <- dplyr::bind_cols(
  dados_fito %>% dplyr::select(parc, spp),
  dap_mat
)

# Rename to dap1..dapN for fitoR auto-detection
colnames(dados_fito)[-(1:2)] <- paste0("dap", seq_len(ncol(dados_fito) - 2))

# Run fitoR (area per plot in m2)
fito_res <- fitoR(dados_fito, area_m2_per_plot = 2500)
fito_res <- data.frame(spp = rownames(fito_res), fito_res, row.names = NULL)

knitr::kable(fito_res, row.names = FALSE) %>%
  kableExtra::kable_styling(full_width = TRUE, position = "center", fixed_thead = TRUE)

# Remove dead/unknown taxa if they appear in spp (optional; adjust to your naming)
dados_fito_clean <- dados_fito %>%
  dplyr::filter(!tolower(spp) %in% c("morta", "indeterminada", "ni"))

fito_res_clean <- fitoR(dados_fito_clean, area_m2_per_plot = 2500)
fito_res_clean <- data.frame(spp = rownames(fito_res_clean), fito_res_clean, row.names = NULL)

# Plot top species by VI components (DR, FR, DoR)
n_sp_plot <- 20

tabela_fito_long <- fito_res_clean %>%
  dplyr::mutate(spp = forcats::fct_reorder(spp, VI)) %>%
  dplyr::select(spp, DR, FR, DoR) %>%
  dplyr::slice_head(n = n_sp_plot) %>%
  tidyr::pivot_longer(cols = c(DR, FR, DoR), names_to = "Parameter", values_to = "Value")

tabela_fito_long %>%
  dplyr::group_by(spp) %>%
  dplyr::mutate(
    Parameter = forcats::fct_relevel(Parameter, "DR", "FR", "DoR"),
    label_y = cumsum(Value)
  ) %>%
  ggplot(aes(x = spp, y = Value, fill = Parameter)) +
  geom_col() +
  geom_text(aes(y = label_y - (Value / 2), label = paste0(round(Value), "%")),
            size = 3, colour = "white") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(
    title = "Species importance value components (top species)",
    subtitle = paste0("Top ", n_sp_plot, " species"),
    x = NULL, y = NULL
  )
