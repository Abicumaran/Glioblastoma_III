setwd("C:/Users/uabic/Desktop")

# Load necessary libraries
library(igraph)
library(dplyr)

# Read the adjacency matrix
adj_matrix <- read.csv("IDHWT_bdm_results_cellrouter.csv", row.names = 1, check.names = FALSE)
adj_matrix <- abs(adj_matrix)
# Convert the adjacency matrix to a graph object
graph <- graph.adjacency(as.matrix(adj_matrix), mode = "undirected", diag = FALSE)

# Compute centrality measures
betweenness_cent <- betweenness(graph)
closeness_cent <- closeness(graph)
eigenvector_cent <- evcent(graph)$vector

# Combine the centrality measures into a data frame
centrality_df <- data.frame(
  Gene = names(betweenness_cent),
  Betweenness = betweenness_cent,
  Closeness = closeness_cent,
  Eigenvector = eigenvector_cent
)

# Select the highest five values for each centrality measure
top_betweenness <- centrality_df %>%
  arrange(desc(Betweenness)) %>%
  slice(1:5)

top_closeness <- centrality_df %>%
  arrange(desc(Closeness)) %>%
  slice(1:5)

top_eigenvector <- centrality_df %>%
  arrange(desc(Eigenvector)) %>%
  slice(1:5)

# Combine the results
top_centrality <- bind_rows(
  top_betweenness %>% mutate(Centrality = "Betweenness"),
  top_closeness %>% mutate(Centrality = "Closeness"),
  top_eigenvector %>% mutate(Centrality = "Eigenvector")
)

# Save the results to a CSV file
write.csv(top_centrality, "Top_Centrality_Measures.csv", row.names = FALSE)

print("Top centrality measures saved to 'Top_Centrality_Measures.csv'")

