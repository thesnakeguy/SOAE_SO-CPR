library(dplyr)
library(cowplot)

# Read data if not already loaded
cpr_data_input <- read.table("SO-CPR_raw_download.txt", sep = "\t")
taxon_df_merged <- read.table("SO-CPR_ID_TaxonAnnotations.txt", sep = "\t")

# Store metadata column names
colnames_metadata <- c("Tow_Number", "Ship_Code", "Time", "Date", 
                       "Month", "Year", "Season", "Latitude", 
                       "Longitude", "Segment_No.", "Segment_Length", 
                       "Total.abundance", "Phytoplankton_Colour_Index", 
                       "Fluorescence", "Salinity", "Water_Temperature", 
                       "Photosynthetically_Active_Radiation")

# Make dataframe to work with
Total_abundance_df <- cpr_data_input %>%
  mutate(
    Total_Plankton_Segment_corr = Total.abundance / Segment_Length,
    Water_Temperature = as.numeric(gsub("-", NA, Water_Temperature))
  )

#############################################
### 1. Total zooplankton abundance index ####
#############################################

effort_df <- Total_abundance_df %>%
  count(Year, name = "n_segments")

year_breaks <- sort(unique(Total_abundance_df$Year))

p_abund <- ggplot(
  Total_abundance_df,
  aes(x = Year, y = Total_Plankton_Segment_corr)
) +
  geom_jitter(
    aes(colour = Water_Temperature),
    width = 0.2,
    alpha = 0.5,
    size = 1
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 10),
    se = FALSE,
    colour = "black"
  ) +
  labs(
    y = "Total plankton (segment corrected)",
    colour = "Water temperature (°C)"
  ) +
  scale_colour_viridis_c(
    option = "plasma",
    na.value = "grey70"
  ) +
  scale_x_continuous(breaks = year_breaks) +
  theme_classic() +
  labs(x = NULL) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    legend.position = "top"
  )

p_effort <- ggplot(
  effort_df,
  aes(x = Year, y = n_segments)
) +
  geom_col(fill = "grey70") +
  labs(y = "Number of segments", x = "Year") +
  scale_x_continuous(breaks = year_breaks) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0.2),
    axis.title.x = element_text(vjust = -0.2)
  )

abundance_plot <- plot_grid(
  p_abund,
  p_effort,
  ncol = 1,
  align = "v",
  rel_heights = c(3, 1)
)
abundance_plot


#######################################
### 2. Salpidae/Euphausiidae index ####
#######################################

Salp.IDs <- taxon_df_merged %>% filter(family == "Salpidae") %>% pull(ID)
Euphausid.IDs <- taxon_df_merged %>% filter(family == "Euphausiidae") %>% pull(ID)
SalpEuph.data <- Total_abundance_df %>%
  dplyr::select(all_of(c(colnames_metadata, Euphausid.IDs, Salp.IDs))) %>%
  mutate(
    Salp.abundance = rowSums(across(all_of(Salp.IDs)), na.rm = TRUE),
    Euph.abundance = rowSums(across(all_of(Euphausid.IDs)), na.rm = TRUE),
    Water_Temperature = as.numeric(Water_Temperature)
  )

year_breaks <- sort(unique(SalpEuph.data$Year))

AnnualBalance <- SalpEuph.data %>% 
  group_by(Year) %>% 
  summarise(
    SalpEuphIndex = mean(Salp.abundance, na.rm = TRUE) - mean(Euph.abundance, na.rm = TRUE)
  )

SalpEuph <- ggplot(SalpEuph.data, aes(x = Year, y = Salp.abundance - Euph.abundance)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50") +
  geom_jitter(
    aes(colour = Salp.abundance - Euph.abundance > 0),
    width = 0.2,
    alpha = 0.5,
    size = 1) +
  geom_point(
    data = AnnualBalance,
    aes(x = Year, y = SalpEuphIndex, colour = "Annual balance"),
    size = 2
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "firebrick",      # segment-level Euphausiid dominance
      "TRUE"  = "steelblue",      # segment-level Salp dominance
      "Annual balance" = "black"  # annual points
    ),
    labels = c(
      "FALSE" = "Euphausiids dominate",
      "TRUE"  = "Salps dominate",
      "Annual balance" = "Annual balance"
    ),
    name = "Dominance"
  ) +
  labs(y = "Salp – Euphausid abundance",
    x = "Year") +
  scale_x_continuous(breaks = year_breaks) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = -45, hjust = 0.2),
    axis.title.x = element_text(vjust = -0.2),
    legend.text = element_text(size = 10)
  ) +
  coord_cartesian(ylim = c(-50, 50))
SalpEuph


###########################
### Save output images ####
###########################
# Abundance plot
ggsave(
  path = "Figures",
  filename = "abundance_timeseries.png",
  plot = abundance_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Salp/Euph plot
ggsave(
  path = "Figures",
  filename = "Salp-Euph_timeseries.png",
  plot = SalpEuph,
  width = 8,
  height = 6,
  dpi = 300
)
  
  
