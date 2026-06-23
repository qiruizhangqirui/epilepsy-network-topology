# =========================================================================
# Script: S0_Demographic.R
# -------------------------------------------------------------------------
# Description:
#   Loads subject demographic information from Subjects.xlsx and generates
#   demographic and clinical distribution plots.
#
#   1. Age distribution for the two centers (JLH & TJU) and combined (All),
#      comparing Epilepsy Patients vs Healthy Controls using a raincloud plot.
#   2. 100% stacked bar charts for epilepsy patients across Cohorts (All, JLH, TJU):
#      - Epilepsy Type (TLE, EXE, GGE, SeLECTS, AE)
#      - FBTCS (Yes vs No)
#      - Seizure Lateralization (Left, Right, Unclear)
#      - Pathology (UHS, BHS, Lesion, Normal)
#
#   Styling and colors are aligned with S11_pySuStaIn_statisitic.m.
#
# Inputs:
#   - data/Subjects.xlsx
#
# Outputs:
#   - outputs/Demographics/age_distribution_JLH.png
#   - outputs/Demographics/age_distribution_TJU.png
#   - outputs/Demographics/age_distribution_All.png
#   - outputs/Demographics/stacked_epilepsy_type.png
#   - outputs/Demographics/stacked_fbtcs.png
#   - outputs/Demographics/stacked_lateralization.png
#   - outputs/Demographics/stacked_pathology.png
#
# Author: Qirui Zhang, Farber Institute for Neuroscience, Department of Neurology, Thomas Jefferson University
# Date: 06/08/2026
# =========================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =========================
# Part 1: Environment & Paths
# =========================
if (dir.exists("data")) {
  project_root <- getwd()
} else {
  project_root <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
}

data_dir <- file.path(project_root, "data")
out_dir  <- file.path(project_root, "outputs", "Demographics")
xlsx_file <- file.path(data_dir, "Subjects.xlsx")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

if (!file.exists(xlsx_file)) {
  stop(sprintf("Input file not found: %s", xlsx_file))
}

message(sprintf("Loading data from: %s", xlsx_file))
df <- read_xlsx(xlsx_file, sheet = 1) %>% as.data.frame()

# =========================
# Part 2: Custom Raincloud Plot Components
# =========================
# Custom GeomFlatViolin definition for ggplot2 raincloud plot
GeomFlatViolin <- ggproto("GeomFlatViolin", GeomViolin,
  setup_data = function(self, data, params) {
    data <- ggplot2:::ggproto_parent(GeomViolin, self)$setup_data(data, params)
    data$xmin <- data$x
    data$xmax <- data$x + data$violinwidth * (data$xmax - data$xmin)
    data
  }
)

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

# =========================
# Part 3: Age Distribution Plot (JLH vs. TJU)
# =========================
message("Generating Age Distribution raincloud plots...")

# Format Group for plotting
df_age <- df %>%
  filter(!is.na(age), !is.na(site)) %>%
  mutate(
    Group_Label = factor(ifelse(Epi == 1, "Epilepsy Patients", "Healthy Controls"),
                         levels = c("Healthy Controls", "Epilepsy Patients"))
  )

# Peach/Orange for Epilepsy Patients (#D99577), Teal/Green for Healthy Controls (#7FB3A2)
colors_age <- c("Healthy Controls" = "#7FB3A2", "Epilepsy Patients" = "#D99577")

# Function to plot age distribution for a single center or combined
plot_center_age <- function(site_code, site_title, out_filename) {
  if (site_code == "All") {
    df_site <- df_age
  } else {
    df_site <- df_age %>% filter(site == site_code)
  }
  
  # Calculate N counts for this site
  n_ctrl <- sum(df_site$Epi == 0)
  n_pat  <- sum(df_site$Epi == 1)
  
  ctrl_lbl <- sprintf("Healthy Controls (N=%d)", n_ctrl)
  pat_lbl  <- sprintf("Epilepsy Patients (N=%d)", n_pat)
  
  # Map labels with N
  df_site <- df_site %>%
    mutate(
      Group_Label_N = factor(
        ifelse(Epi == 1, pat_lbl, ctrl_lbl),
        levels = c(ctrl_lbl, pat_lbl)
      )
    )
  
  p <- ggplot(df_site, aes(x = Group_Label_N, y = age, fill = Group_Label)) +
    # Half violin (flat violin) nudged to the right (faces upwards in coord_flip)
    geom_flat_violin(position = position_nudge(x = 0.15, y = 0), adjust = 1.2, 
                     color = "black", size = 0.8, alpha = 0.8) +
    # Raw data points nudged to the left with jitter
    geom_point(aes(x = as.numeric(Group_Label_N) - 0.15, color = Group_Label),
               position = position_jitter(width = 0.05, height = 0), 
               size = 1.5, alpha = 0.5) +
    # Narrow boxplot centered
    geom_boxplot(width = 0.1, position = position_nudge(x = 0, y = 0), 
                 outlier.shape = NA, color = "black", size = 0.8, fill = "white", alpha = 0.9) +
    # Custom colors
    scale_fill_manual(values = colors_age) +
    scale_color_manual(values = colors_age) +
    # Horizontal orientation
    coord_flip() +
    # Labels and themes
    labs(title = sprintf("Age Distribution (%s)", site_title), x = "", y = "Years") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1.2, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text.y = element_text(face = "bold", size = 13),
      axis.text.x = element_text(size = 13),
      legend.position = "none" # Hide legend
    )
  
  ggsave(file.path(out_dir, out_filename), plot = p, width = 8, height = 4, dpi = 300)
}

# Generate plots for the two centers and the combined cohort
plot_center_age("JLH", "JLH Center", "age_distribution_JLH.png")
plot_center_age("TJU", "TJU Center", "age_distribution_TJU.png")
plot_center_age("All", "All Centers", "age_distribution_All.png")

message("Saved age_distribution_JLH.png, age_distribution_TJU.png, and age_distribution_All.png")

# =========================
# Part 4: Helper for Stacked Bar Charts
# =========================
# Replicates style of local_stacked_bar.m:
# - Y-axis scaled 0-100%
# - Legend on the right side
# - Thick axes, no gridlines
# - Specific color palettes
plot_stacked_bar <- function(data_df, var_name, fill_colors, legend_labels, levels_order) {
  # Enforce factor levels and ordering
  data_df$Category <- factor(data_df$Category, levels = levels_order)
  data_df$Cohort <- factor(data_df$Cohort, levels = c("All", "JLH", "TJU"))
  
  # Calculate percentages
  stacked_pct <- data_df %>%
    group_by(Cohort, Category) %>%
    summarise(count = n(), .groups = 'drop_last') %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  # Base plot
  p <- ggplot(stacked_pct, aes(x = Cohort, y = percentage, fill = Category)) +
    geom_bar(stat = "identity", width = 0.5, color = NA) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%")) +
    scale_fill_manual(values = fill_colors, labels = legend_labels) +
    labs(title = var_name, x = "", y = "") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.ticks = element_line(size = 1.2, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text = element_text(face = "bold", size = 13),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "right"
    )
  
  return(p)
}

# Patients subset for stacked bar charts
patients <- df %>% filter(Epi == 1)

# Duplicate data to represent cohorts: All, JLH, TJU
prepare_cohort_data <- function(subset_df) {
  d_all <- subset_df %>% mutate(Cohort = "All")
  d_jlh <- subset_df %>% filter(site == "JLH") %>% mutate(Cohort = "JLH")
  d_tju <- subset_df %>% filter(site == "TJU") %>% mutate(Cohort = "TJU")
  bind_rows(d_all, d_jlh, d_tju)
}

# =========================
# Part 5: Stacked Bar Chart 1 - Epilepsy Type
# =========================
message("Generating Epilepsy Type stacked bar chart...")
# Proportions of TLE, EXE, GGE, SeLECTS, AE
# Palette matching S11 and extensions:
# TLE: deep blue (#2659B2), EXE: light blue (#8CB3E6), GGE: teal (#3B9AB2), SeLECTS: sky blue (#78B7C5), AE: purple-blue (#7E60B2)
colors_epi_type <- c("TLE" = "#2659B2", "EXE" = "#8CB3E6", "GGE" = "#3B9AB2", "SeLECTS" = "#78B7C5", "AE" = "#7E60B2")
epi_type_order  <- c("TLE", "EXE", "GGE", "SeLECTS", "AE")

df_epi_type <- prepare_cohort_data(patients) %>%
  rename(Category = Group) %>%
  filter(Category %in% epi_type_order)

p_epi_type <- plot_stacked_bar(
  data_df = df_epi_type,
  var_name = "Epilepsy Type",
  fill_colors = colors_epi_type,
  legend_labels = names(colors_epi_type),
  levels_order = epi_type_order
)

ggsave(file.path(out_dir, "stacked_epilepsy_type.png"), plot = p_epi_type, width = 6, height = 4.5, dpi = 300)

# =========================
# Part 6: Stacked Bar Chart 2 - FBTCS
# =========================
message("Generating FBTCS stacked bar chart...")
# Proportions of Yes vs No (1=Yes, 0=No in Subjects.xlsx)
# Palette matching S11: No: deep red (#B32630), Yes: light red (#F38F8F)
colors_fbtcs <- c("No" = "#B32630", "Yes" = "#F38F8F")
fbtcs_order  <- c("No", "Yes")

# Filter out NA values for FBTCS
df_fbtcs <- prepare_cohort_data(patients) %>%
  filter(!is.na(`FTBTC(Y/N)`)) %>%
  mutate(Category = ifelse(`FTBTC(Y/N)` == 1, "Yes", "No"))

p_fbtcs <- plot_stacked_bar(
  data_df = df_fbtcs,
  var_name = "FBTCS",
  fill_colors = colors_fbtcs,
  legend_labels = names(colors_fbtcs),
  levels_order = fbtcs_order
)

ggsave(file.path(out_dir, "stacked_fbtcs.png"), plot = p_fbtcs, width = 6, height = 4.5, dpi = 300)

# =========================
# Part 7: Stacked Bar Chart 3 - Seizure Lateralization
# =========================
message("Generating Seizure Lateralization stacked bar chart...")
# Proportions of Left, Right, Unclear
# Palette matching S11: Left: deep orange (#D9661A), Right: light orange (#FAA659), Unclear: gray (#B3B3B3)
colors_lat <- c("Left" = "#D9661A", "Right" = "#FAA659", "Unclear" = "#B3B3B3")
lat_order  <- c("Left", "Right", "Unclear")

df_lat <- prepare_cohort_data(patients) %>%
  mutate(Category = case_when(
    lateralization_L == 1 ~ "Left",
    lateralization_R == 1 ~ "Right",
    lateralization_UC == 1 ~ "Unclear",
    TRUE ~ "Unclear"
  ))

p_lat <- plot_stacked_bar(
  data_df = df_lat,
  var_name = "Seizure Lateralization",
  fill_colors = colors_lat,
  legend_labels = names(colors_lat),
  levels_order = lat_order
)

ggsave(file.path(out_dir, "stacked_lateralization.png"), plot = p_lat, width = 6, height = 4.5, dpi = 300)

# =========================
# Part 8: Stacked Bar Chart 4 - Pathology
# =========================
message("Generating Pathology stacked bar chart...")
# Proportions of UHS, BHS, Lesion, Normal
# Palette matching S11: UHS: deep purple (#733399), BHS: light purple (#B373CC), Lesion: light violet (#E0C0ED), Normal: gray (#BFBFBF)
colors_path <- c("UHS" = "#733399", "BHS" = "#B373CC", "Lesion" = "#E0C0ED", "Normal" = "#BFBFBF")
path_order  <- c("UHS", "BHS", "Lesion", "Normal")

df_path <- prepare_cohort_data(patients) %>%
  mutate(Category = case_when(
    Pathology_UHS == 1 ~ "UHS",
    Pathology_BHS == 1 ~ "BHS",
    Pathology_Lesion == 1 ~ "Lesion",
    Pathology_Normal == 1 ~ "Normal",
    TRUE ~ "Normal"
  ))

p_path <- plot_stacked_bar(
  data_df = df_path,
  var_name = "Pathology",
  fill_colors = colors_path,
  legend_labels = names(colors_path),
  levels_order = path_order
)

ggsave(file.path(out_dir, "stacked_pathology.png"), plot = p_path, width = 6, height = 4.5, dpi = 300)

message("Demographic and clinical visualization script execution completed.")
message(sprintf("All figures saved to: %s", normalizePath(out_dir)))
