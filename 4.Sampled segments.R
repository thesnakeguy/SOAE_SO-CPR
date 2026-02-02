### DESCRIPTION #######################################################################
## This code can be used to visualize the segments that were sampled in a specific year
#######################################################################################


# Load required libraries
library(sf)
library(rnaturalearth)
library(ggplot2)
library(SOmap)

# Parameters
TargetYear = "2022"

# Read data 
cpr_data_input2 <- read.table("SO-CPR_raw_download.txt", sep = "\t")
comm_df <- read.table("Community_assignments.txt", sep = "\t", header = TRUE, comment.char = "")

# Make unique segment identifier
cpr_data_input2 <- cpr_data_input2 %>% mutate(Segment = as.character(paste(Ship_Code, Tow_Number, Date, Time, sep = "_")))

# Prepare background data
data_background <- cpr_data_input2 %>%
  mutate(
    Longitude = as.numeric(Longitude),
    Latitude = as.numeric(Latitude)
  ) %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  sf::st_transform(crs = 3031)

# Prepare target year data
data_year <- data_background %>% filter(Year == TargetYear)

# Get land data
land <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 3031)

# Get SOmap and fronts data
fronts <- SOmap_data$fronts_orsi %>% st_as_sf %>% st_transform(crs = 3031)

################################
### PLOTTING ###################
################################

# Plot sampled segments with community assignment information
# Join community assignment info from 3.Network_analysis with the working data
comm_df$Segment <- as.character(comm_df$Segment)
comm_df$Comm_colors <- as.factor(comm_df$Comm_colors)

# Merge community data with year data
data_year_comm <- data_year %>%
  left_join(comm_df, by = c("Segment"))
data_year_comm$Community <- as.factor(data_year_comm$Community)
comm_color_map <- setNames(comm_df$Comm_colors, comm_df$Community)

# Create a named vector of colors
comm_color_map <- setNames(comm_df$Comm_colors, comm_df$Community)
sampled_segments_comm <- ggplot() +
  geom_sf(data = land, fill = "antiquewhite", color = "darkgrey", size = 0.2) +
  geom_sf(data = data_background, color = "grey80", size = 0.5, alpha = 0.2) +
  geom_sf(data = data_year_comm, aes(color = Community), size = 3) +
  geom_sf(data = fronts, color = "steelblue", size = 1, pch = 2) +
  scale_color_manual(values = comm_color_map, name = "Community") +
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
    title = paste0("Plankton community network (", TargetYear, ")"),
    subtitle = paste0(length(data_year$Segment), " Segments sampled in ", length(unique(data_year$Tow_Number)), " tows"),
    x = "Longitude (EPSG:3031)",
    y = "Latitude (EPSG:3031)"
  )
sampled_segments_comm

####################
### Save output ####
####################

ggsave(
  path = "Figures",
  filename = paste0("Sampled_Segments_comm_", TargetYear, ".png"),
  plot = sampled_segments,
  width = 8,
  height = 6,
  dpi = 300
)
