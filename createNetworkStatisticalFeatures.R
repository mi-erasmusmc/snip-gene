# Calculate the network statistic features

# Load the required packages
require(igraph)
require(orca)

if(config$recalculateNetwork){

  # Prepareer data voor count functie orca
  nodes = data.frame(Gene = unique(c(rels$sn, rels$on)))
  nodes$ID = 1:nrow(nodes)
  
  ppi_undirected = graph_from_edgelist(as.matrix(rels[,c("sn", "on")]), directed = F)
  ppi_undirected = simplify(ppi_undirected)
  undirected_edges = as.data.frame(as_edgelist(ppi_undirected), stringsAsFactors = F)
  undirected_edges = merge(undirected_edges, nodes, by.x = "V2", by.y = "Gene")
  undirected_edges = merge(undirected_edges, nodes, by.x = "V1", by.y = "Gene")
  
  # Identificeer graphlets conform Hocevar en Demsar
  if(!config$long){
    graphlets = as.data.frame(count4(undirected_edges[,c("ID.x", "ID.y")]))
  } else {
    graphlets = as.data.frame(count5(undirected_edges[,c("ID.x", "ID.y")]))
  }
  colnames(graphlets) = paste0("Graphlet_", colnames(graphlets))
  rownames(graphlets) = nodes$Gene
  
  # Other metrics
  eigen_values = as.data.frame(eigen_centrality(ppi_undirected)$vector)
  betweenness_values = as.data.frame(betweenness(ppi_undirected, directed = F))
  
  other = merge(eigen_values, betweenness_values, by = "row.names", all = T)
  colnames(other) = c("Gene", "Eigen_centrality", "Betweenness")
  
  network_features = merge(other, graphlets, by.x = "Gene", by.y = "row.names", all = T)
  
  # Large scale cluster identification
  ppi_clusters = cluster_louvain(simplify(ppi_undirected))
  ppi_clusters = data.frame(Gene = ppi_clusters$names, Louvain_cluster = ppi_clusters$membership)
  ppi_clusters$Louvain_cluster = paste0("Cluster_", ppi_clusters$Louvain_cluster)
  
  network_features = merge(network_features, ppi_clusters, by = "Gene", all = T)
  
  if(config$long){
    betweenness_clusters = cluster_edge_betweenness(ppi_undirected) # Duurt ook lang, maar is wel relevant (zie publicatie)
    ppi_clusters2 = data.frame(Gene = as.vector(V(g)), Betweenness_cluster = betweenness_clusters$membership)
    ppi_clusters2$Betweenness_cluster = paste0("Cluster_", ppi_clusters2$Betweenness_cluster)
    network_features = merge(network_features, ppi_clusters2, by = "Gene", all = T)
  }
  write.csv2(network_features, paste0("Raw data files/Network metric features calculated on ", todays_date, ".csv"), row.names = F)
} else {
  network_features = read.csv2(paste0("Raw data files/Network metric features calculated on ", config$Data.date, ".csv"), stringsAsFactors = F)
}