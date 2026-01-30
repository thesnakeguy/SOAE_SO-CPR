

# Load the CPR data from 0.Download_SO-CPR_data.R
cpr_data_input <- read.table("SO-CPR_raw_download.txt", sep = "\t")
# Column names of metadata and plankton species
colnames_metadata <- c("Tow_Number", "Ship_Code", "Time", "Date", 
                       "Month", "Year", "Season", "Latitude", 
                       "Longitude", "Segment_No.", "Segment_Length", 
                       "Total.abundance", "Phytoplankton_Colour_Index", 
                       "Fluorescence", "Salinity", "Water_Temperature", 
                       "Photosynthetically_Active_Radiation")
colnames_species <- setdiff(colnames(cpr_data_input), colnames_metadata)


#################################
# Function to parse taxon ID ####
#################################

life_stages <- c("egg", "nauplius", "metanauplius", "zoea", "megalopa", "phyllosoma",
                 "juv", "natant", "larvae", "small", "calyptopis", "furcilia", "cyprid",
                 "nectophore")
`%ni%` <- Negate(`%in%`)

parse_taxa_id <- function(id) {
  parts <- str_split(id, "\\.")[[1]]
  
  # Handle Egg cases
  if (str_detect(id, "^Egg\\.?")) {
    return(list(ID = id, LifeStage = "egg"))
  }
  
  # Handle life stage codes (C1, F2, etc. eg: Euphausia.superba.F3)
  if (str_detect(id, "\\.[A-Za-z]+\\.[CF][0-9]$")) {
    parts <- str_split(id, "\\.", n = 3)[[1]]
    return(list(ID = id, Genus = parts[1], Species = parts[2], LifeStage = parts[3]))
  }
  
  # Handle species with life stage (eg: Euphausia.triacantha.calyptopis)
  if (str_detect(id, "^[A-Z][a-z]+\\.[a-z]+\\.[a-z]+$") & parts[3] %in% life_stages) {
    return(list(ID = id, Genus = parts[1], Species = parts[2], LifeStage = parts[3]))
  }
  
  # Case: "Genus..Subgenus..species" -> Keep Genus and Species
  if (str_detect(id, "^[A-Z][a-z]+\\.\\.[A-Z][a-z]+\\.\\.[a-z]+$")) {
    parts <- str_split(id, "\\.\\.")[[1]]
    return(list(ID = id, Genus = parts[1], Species = parts[3]))
  }
  
  # Handle genus level cases (Genus.sp)
  if (str_detect(id, "\\.sp\\.$")) {
    genus <- str_remove(id, "\\.sp\\.?$")
    return(list(ID = id, Genus = genus, Qualifier = "sp"))
  }
  
  # Handle genus with life stage (eg: Thysanoessa.sp..furcilia)
  if (str_detect(id, "\\.sp\\.\\.[a-z]+$") & (tolower(parts[length(parts)]) %in% life_stages)) {
    parts <- str_split(id, "\\.sp\\.\\.", n = 2)[[1]]
    return(list(ID = id, Genus = parts[1], LifeStage = parts[2], Qualifier = "sp"))
  }
  
  # Handle higher taxon with life stage (eg: bryozoa.larvae, Decapoda.megalopa)
  if (tolower(parts[length(parts)]) %in% life_stages & length(parts == 2)) {
    return(list(ID = id, HigherTaxon = parts[1], LifeStage = parts[2]))
  }
  
  # Handle higher taxon with life stages and "indet" (eg: Calanoida.indet..small. or Copepoda.nauplius.indet)
  if (str_detect(id, "\\.indet\\.\\.[a-z]")) {
    parts <- str_split(id, "\\.indet\\.\\.", n = 2)[[1]]
    return(list(ID = id, HigherTaxon = parts[1], LifeStage = parts[2], Qualifier = "indet"))
  }
  if (str_detect(id, "[A-Z][a-z]+\\.[a-z]+\\.indet$") & parts[2] %in% life_stages) {
    return(list(ID = id, HigherTaxon = parts[1], LifeStage = parts[2], Qualifier = "indet"))
  }
  
  # Case: "Taxon.indet" -> HigherTaxon, Qualifier
  if (str_detect(id, "\\.indet$")) {
    taxon <- str_remove(id, "\\.indet$")
    return(list(ID = id, HigherTaxon = taxon, Qualifier = "indet"))
  }
  
  # Handle subspecies (eg: Clione.limacina.antarctica)
  if (str_detect(id, "^[A-Z][a-z]+\\.[a-z]+\\.[a-z]+$") & parts[3] %ni% life_stages) {
    return(list(ID = id, Genus = parts[1], Species = parts[2]))
  }
  
  # Handle species with "var" (eg: Euphausia.similis.var.armata)
  if (length(grep("var", id)) > 0) {
    return(list(ID = id, Genus = parts[1], Species = parts[2]))
  }
  
  # Default genus-species case
  if (str_detect(id, "^[A-Z][a-z]+\\.[a-z]+$") & (tolower(parts[length(parts)]) %ni%  life_stages)) {
    parts <- str_split(id, "\\.")[[1]]
    return(list(ID = id, Genus = parts[1], Species = parts[2]))
  }
  
  # Default case
  return(list(ID = id, RawID = id))
}




#################################################
# Function to retrieve hierarchical taxonomy ####
#################################################


process_taxon <- function(taxon) {
  
  # Initialize taxon hierarchy with NA for all ranks
  columns <- c("ID", "genus", "species", "family", "suborder", "order", "subclass", "class", "subphylum",
               "phylum", "infrakingdom", "subkingdom", "kingdom",
               "LifeStage", "Qualifier")
  taxon_hierarchy <- setNames(rep(NA, length(columns)), columns)
  
  
  
  # Retrieve information based on ID (Genus + Species, Genus, HigherTaxon, LifeStage)
  AphiaID <- NA  # Initialize AphiaID to NA in case all conditions fail
  tryCatch({
    if (!is.null(taxon$Genus) && !is.null(taxon$Species)) {
      # Genus + Species: Query based on both genus and species
      AphiaID <- wm_name2id(paste(taxon$Genus, taxon$Species))
    } else if (!is.null(taxon$Genus)) {
      # Genus only: Query based on genus
      AphiaID <- wm_name2id(taxon$Genus)
    } else if (!is.null(taxon$HigherTaxon)) {
      # HigherTaxon only: Query based on higher taxon (use it as a higher taxonomic rank)
      AphiaID <- wm_name2id(taxon$HigherTaxon)
    }
  }, error = function(e) {
    message("Error in wm_name2id for taxon: ", taxon$Genus, " ", taxon$Species, " - ", e$message)
    AphiaID <<- NA  # Assign NA if error occurs
  })
  
  # If AphiaID is still NA, return early with taxon_hierarchy
  if (is.na(AphiaID)) {
    return(taxon_hierarchy)
  }
  
  # Fetch classification if AphiaID is valid
  if (!is.null(AphiaID) && AphiaID != "") {
    classification <- wm_classification(AphiaID)
    
    # Loop through the ranks and fill taxon hierarchy dynamically
    for (i in seq_along(classification[[2]])) {
      rank <- tolower(classification[[2]][i])  
      name <- classification[[3]][i]
      
      if (rank %in% names(taxon_hierarchy)) {
        taxon_hierarchy[[rank]] <- name
      }
    }
  }
  
  # Handle life stage if available and add the original ID
  taxon_hierarchy$ID <- taxon$ID
  if (!is.null(taxon$LifeStage)) {
    taxon_hierarchy$LifeStage <- taxon$LifeStage
  }
  if (!is.null(taxon$Qualifier)) {
    taxon_hierarchy$Qualifier <- taxon$Qualifier  # Life stage could be in the qualifier
  }
  
  return(taxon_hierarchy)
}

# List of taxonomic IDs
taxon_ids <- colnames_species #these are the ID's of the SO-CPR raw data that represent species/Lifestage mix

# Parse IDs
parsed_taxa <- map(taxon_ids, parse_taxa_id)

# Apply function to each taxon in parsed_taxa and bind results into a data frame
taxon_df <- bind_rows(lapply(parsed_taxa, process_taxon))

# Convert to data frame for easy viewing
taxon_df <- as.data.frame(taxon_df)

# Add rows of IDs that are horribly picked to automate the process
taxon_extra <- data.frame(
  ID = c("Appendicularia.indet", "Clione.sp.", "Ctenophora.indet", "Hyperia.sp.", "Spongiobranchaea.australis", "Themisto.sp."),
  genus = c("Appendicularia", "Clione", NA, "Hyperia", "Spongiobranchaea", "Themisto"),
  species = c(NA, NA, NA, NA, "australis", NA),
  family = c(NA, "Clionidae", NA, "Hyperiidae", "Pneumodermatidae", "Hyperiidae"),
  suborder = c(NA, NA, NA, NA, NA, NA),
  order = c(NA, "Pteropoda", NA, "Amphipoda", "Pteropoda", "Amphipoda"),
  subclass = c(NA, NA, NA, NA, NA, NA),
  class = c("Appendicularia", "Gastropoda", NA, "Malacostraca", "Gastropoda", "Malacostraca"),
  subphylum = c(NA, NA, NA, NA, NA, NA),
  phylum = c("Chordata", "Mollusca", "Ctenophora", "Arthropoda", "Mollusca", "Arthropoda"),
  infrakingdom = c(NA, NA, NA, NA, NA, NA),
  subkingdom = c(NA, NA, NA, NA, NA, NA),
  kingdom = c("Animalia", "Animalia", "Animalia", "Animalia", "Animalia", "Animalia"),
  LifeStage = c(NA, NA, NA, NA, NA, NA),
  Qualifier = c("indet", "sp", "indet", "sp", NA, "sp"),
  stringsAsFactors = FALSE
) 

# Merge with your existing dataframe
names(taxon_extra) <- names(taxon_df)
taxon_df_merged <- rbind(taxon_df, taxon_extra)
taxon_df_merged <- taxon_df_merged[rowSums(!is.na(taxon_df_merged)) > 0, ]

# Print final data frame
view(taxon_df_merged, n = Inf)

# Save the dataframe
write.table(taxon_df_merged, "SO-CPR_ID_TaxonAnnotations.txt", quote = FALSE, sep = "\t")
