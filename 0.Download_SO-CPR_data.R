remotes::install_github("AustralianAntarcticDivision/blueant")
library(blueant)

## Read CPR plankton data using BlueAnt ####
# Configure data folder
dir.create("SO-CPR_RawData_BlueAnt/")
cpr_plankton_folder <- "SO-CPR_RawData_BlueAnt/"
cpr_plankton <- bb_config(local_file_root = cpr_plankton_folder)

# Select and add CPR plankton dataset from BlueAnt
SO_plankton <- sources("Southern Ocean Continuous Plankton Recorder")
cpr_plankton <- cpr_plankton %>% bb_add(SO_plankton)

# Synchronize/download the dataset
status <- bb_sync(cpr_plankton)

# Load the CPR data
myfiles <- status$files[[1]]
cpr_data_input <- read.csv(myfiles$file[grepl("AADC", myfiles$file)])

# Column names of metadata and plankton species
colnames_metadata <- c("Tow_Number", "Ship_Code", "Time", "Date", 
                       "Month", "Year", "Season", "Latitude", 
                       "Longitude", "Segment_No.", "Segment_Length", 
                       "Total.abundance", "Phytoplankton_Colour_Index", 
                       "Fluorescence", "Salinity", "Water_Temperature", 
                       "Photosynthetically_Active_Radiation")
colnames_species <- setdiff(colnames(cpr_data_input), colnames_metadata)

# Save the dataset
write.table(cpr_data_input, file = "SO-CPR_raw_download.txt", sep="\t", quote = FALSE)
