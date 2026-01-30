# Load required libraries
library(sf)
library(rnaturalearth)
library(ggplot2)

# Parameters
TargetYear = "2022"

# Read data if not already loaded
cpr_data_input <- read.table("SO-CPR_raw_download.txt", sep = "\t")

# Prepare background data
data_background <- cpr_data_input %>%
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude = as.numeric(Latitude)
  ) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = 3031)

# Prepare target year data
data_year <- data_background %>% filter(Year == TargetYear)

# Get land data
land <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 3031)

# Get SOmap and fronts data
fronts <- SOmap_data$fronts_orsi %>% st_as_sf %>% st_transform(crs = 3031)

# Plot distribution
sampled_segments <- ggplot() +
  geom_sf(data = land, fill = "antiquewhite", color = "darkgrey", size = 0.2) +
  geom_sf(data = data_background, color = "grey80", size = 0.5, alpha = 0.2) +
  geom_sf(data = data_year, color = "firebrick", size = 1) +
  geom_sf(data = fronts, color = "steelblue", size = 1, pch = 2) +
  coord_sf(crs = st_crs(3031), ylim = c(-6500000, 6500000), xlim = c(-6500000, 6500000)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "aliceblue"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "gray90"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(
    title = paste0("Segments sampled in ",TargetYear),
    subtitle = "Spatial distribution with oceanic fronts",
    x = "Longitude (EPSG:3031)",
    y = "Latitude (EPSG:3031)"
  )
sampled_segments

# Save output
ggsave(
  path = "Figures",
  filename = paste0("Sampled_Segments_", TargetYear, ".png"),
  plot = sampled_segments,
  width = 8,
  height = 6,
  dpi = 300
)
