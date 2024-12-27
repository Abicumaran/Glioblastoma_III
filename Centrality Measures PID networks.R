
# Load necessary libraries
library(igraph)
library(data.table)

# Function to calculate centrality measures and save to CSV
calculate_and_save_centralities <- function(file_path, output_dir) {
  # Read the data
  data <- fread(file_path, header = FALSE)
  colnames(data) <- c("Gene1", "Gene2", "Weight")
  
  # Create the graph
  graph <- graph_from_data_frame(data, directed = FALSE)
  
  # Calculate centrality measures
  eigen_centrality <- eigen_centrality(graph, directed = FALSE)$vector
  betweenness_centrality <- betweenness(graph, directed = FALSE)
  closeness_centrality <- closeness(graph, mode = "all")
  
  # Create a data frame with the results
  centrality_measures <- data.frame(
    Gene = names(eigen_centrality),
    Eigenvector_Centrality = eigen_centrality,
    Betweenness_Centrality = betweenness_centrality,
    Closeness_Centrality = closeness_centrality
  )
  
  # Define the output file name
  output_file <- paste0(output_dir, "/", basename(file_path), "_centralities.csv")
  
  # Save the results to CSV
  fwrite(centrality_measures, file = output_file, row.names = FALSE)
}

# List of input files and output directory
input_files <- c("PID_selected_HVGIDHWT.txt", "PID_selected_HVGK27M.txt", 
                 "PID_selected_scEpathIDHWT.txt", "PID_selected_scEpathK27M.txt", 
                 "PID_selected_TopFeatIDHWT.txt", "PID_selected_TopFeatK27M.txt")

output_directory <- "C:/Users/uabic/Desktop/CentralityResults"

# Ensure output directory exists
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Process each file
for (file in input_files) {
  calculate_and_save_centralities(file, output_directory)
}
