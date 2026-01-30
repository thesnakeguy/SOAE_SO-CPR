### DESCRIPTION ############################################################################
# This code can be used to explore similarities in plankton communities between 
# standardized transects (segments) in the SO-CPR dataset. It applies a network analysis 
# on Bray-Curtis distances, relying on the Igraph package for visualization 
############################################################################################

# Load required libraries
library(igraph)
library(dplyr)
library(tidyr)
library(vegan)
library(backbone)


# Read data if not already loaded
cpr_data_input <- read.table("SO-CPR_raw_download.txt", sep = "\t")
taxon_df_merged <- read.table("SO-CPR_ID_TaxonAnnotations.txt", sep = "\t")

# Make unique segment identifier
cpr_data_input <- cpr_data_input %>% mutate(Segment_uniq = paste(Ship_Code, Tow_Number, Date, Time, sep = "_"))

# Select only the taxa columns and build segments dataset
colnames_metadata <- c("Tow_Number", "Ship_Code", "Time", "Date", 
                       "Month", "Year", "Season", "Latitude", 
                       "Longitude", "Segment_No.", "Segment_Length", 
                       "Total.abundance", "Phytoplankton_Colour_Index", 
                       "Fluorescence", "Salinity", "Water_Temperature", 
                       "Photosynthetically_Active_Radiation", "Segment_uniq")

taxa_data <- cpr_data_input[, -match(colnames_metadata, names(cpr_data_input))]
segments <- cbind(cpr_data_input$Segment_uniq, cpr_data_input$Latitude, cpr_data_input$Year, taxa_data)
colnames(segments)[1:3] <- c("Segment", "Latitude", "Year")
full_data <- cbind(segments, taxa_data)

# Select only taxa columns for similarity calculation
full_data[is.na(full_data)] <- 0

# Subset data if needed for testing purposes
sub_data <- full_data[1:10000,]

# Compute similarity matrix (using Bray-Curtis as an example)
distance_matrix <- vegdist(sub_data[,-c(1:3)], method = "bray")  # Dissimilarity matrix
distance_matrix <- as.matrix(distance_matrix)
similarity_matrix <- 1 - distance_matrix  # Convert to similarity

net <- graph_from_adjacency_matrix(
  similarity_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

V(net)$name <- sub_data$Segment 
V(net)$year <- sub_data$Year
V(net)$latitude <- sub_data$Latitude 


# Network backbone extraction
backbone_net <- backbone_from_weighted(net, model = "disparity", alpha = 0.05)

# Remove vertices that have less than X edges
deg_backbone_net <- degree(backbone_net)
keep <- which(deg_backbone_net >= 50)
backbone_net_pruned <- induced_subgraph(backbone_net, keep)

# Plot net 
plot(backbone_net_pruned,
     vertex.size = 2,
     vertex.label = NA,
     vertex.label.cex = 0.5,
     vertex.color = rainbow(length(unique(V(net)$year)), alpha = 0.2)[as.factor(V(net)$year)],
     edge.width = 0.2, #alternatively: E(net)$weight * 5
     edge.color = "grey80",
     layout = layout_with_fr
)

# Get unique latitudes and assign colors
color_names <- unique(V(net)$year)
colors <- rainbow(length(color_names))
names(colors) <- color_names

# Add legend
legend("topright",
       legend = round(color_names, 2), 
       fill = colors,
       title = "Year",
       bty = "n",  
       cex = 1   
)

