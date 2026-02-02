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
library(RColorBrewer)
library(randomcoloR)

###################
### Parameters ####
###################

TargetYear = "2022"

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
segments <- cbind(cpr_data_input$Segment_uniq,
                  cpr_data_input$Tow_Number,
                  cpr_data_input$Latitude,
                  cpr_data_input$Longitude,
                  cpr_data_input$Year,
                  taxa_data)
colnames(segments)[1:5] <- c("Segment", "Tow_Number", "Latitude", "Longitude", "Year")

# Filter data for TargetYear
sub_data <- segments %>% filter(Year == TargetYear)
dim(sub_data)

# Select only taxa columns for similarity calculation
sub_data[is.na(sub_data)] <- 0


# Compute similarity matrix (using Bray-Curtis as an example)
distance_matrix <- vegdist(sub_data[,-c(1:5)], method = "bray")  # Dissimilarity matrix
distance_matrix <- as.matrix(distance_matrix)
similarity_matrix <- 1 - distance_matrix  # Convert to similarity

net <- graph_from_adjacency_matrix(
  similarity_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)
V(net)$Segment <- sub_data$Segment 
V(net)$Tow_Number <- sub_data$Tow_Number
V(net)$Latitude <- sub_data$Latitude
V(net)$Longitude <- sub_data$Longitude

# Network backbone extraction
backbone_net <- backbone_from_weighted(net, model = "disparity", alpha = 0.05)
V(backbone_net)$Segment <- V(net)$Segment
V(backbone_net)$Tow_Number <- V(net)$Tow_Number
V(backbone_net)$Latitude <- V(net)$Latitude
V(backbone_net)$Longitude <- V(net)$Longitude

# Keep only the largest connected component
components <- components(backbone_net)
largest_component <- which.max(components$csize)
keep <- which(components$membership == largest_component)
backbone_net_pruned <- induced_subgraph(backbone_net, keep)
deg_backbone_net_pruned <- degree(backbone_net_pruned)
backbone_net_pruned <- induced_subgraph(backbone_net_pruned, which(deg_backbone_net_pruned >= 5))



##############################
### Plotting #################
##############################
### A. Plot with tow number colors

Net_selected <- net # one of net, backbone_net, backbone_net_pruned

# Set colors with Tow_Number
tows <- sort(unique(V(Net_selected)$Tow_Number))
n_tows <- length(tows)
tow_cols <- setNames(
  distinctColorPalette(n_tows)[seq_len(n_tows)],
  tows
)

# Define graph layout
fr_layout <- layout_with_fr(Net_selected)

# Plot net with tow number as color code
plot(Net_selected,
     vertex.size = 5,
     vertex.label = NA,
     vertex.label.cex = 0.5,
     vertex.color = tow_cols[as.character(V(Net_selected)$Tow_Number)],
     edge.width = 0.2, 
     edge.color = "grey80",
     main = NULL,
     sub = paste0(length(unique(V(Net_selected)$Tow_Number)), " Tows visualized"),
     layout = fr_layout
)

# Add legend
legend(
  "topright",
  legend = names(tow_cols),
  col = tow_cols,
  pch = 16,
  pt.cex = 1.5,
  bty = "n",
  title = "Tow number"
)


### B. Plot net with community detection
# Community detection algorithm
comm <- cluster_walktrap(Net_selected)
# Create dataframe with assignment and colors (to be reused in 4.Sampled segments.R)
comm_df <- data.frame(
  Segment = as.character(V(Net_selected)$Segment),  
  Community = as.factor(membership(comm))
)
n_comm <- length(unique(comm_df$Community))
comm_cols <- setNames(
  distinctColorPalette(n_comm),
  as.character(1:n_comm)  
)
comm_df$Comm_colors <- comm_cols[as.character(membership(comm))]

plot(Net_selected,
     vertex.size = 5,
     vertex.label = NA,
     vertex.label.cex = 0.5,
     vertex.color = comm_df$Comm_colors,
     edge.width = 0.2,
     edge.color = "grey80",
     main = NULL,
     sub = paste("Walktrap algorithm -", length(comm), "communities detected"),
     layout = fr_layout
)


########################
### Save the output ####
########################

write.table(comm_df, "Community_assignments.txt", quote = FALSE, sep = "\t")
png(
  filename = "Figures/Community_Network_Yearcoded.png",
  width = 800,     
  height = 600,    
  units = "px",    
  res = 300        
)
# Now plot the image
dev.off()
