# Load the required packages
require(mongolite)
require(data.table)
require(httr)
require(rjson)
require(igraph)
require(tidyr)

# Parameters
fold_change_cutoff = 0.1
FDR_cutoff = 0.001
proteinOnly = T

# Select the relevant genes
ngs = read.csv("GeneExpressionData/Galaxy37-[edgeR_DGE_on_2__design_matrix_prostate_unpaired.txt_-_differentially_expressed_genes].tabular.annotated.txt", 
                           stringsAsFactors = F, sep = "\t", row.names = 1)
colnames(ngs) = c("Ensembl_ID", "Name", "log2FC", "logCPM", "LR", "p_value", "FDR")
ngs$Ensembl_ID_no_Dot = sapply(ngs$Ensembl_ID, function(x){return(unlist(strsplit(x, "\\."))[1])})
ngs$FC = 2^ngs$log2FC

# Get the locations of the genes
gene_data = read.csv("GeneExpressionData/ENSEMBL gene locations downloaded from BioMart on 13-05-2019.txt", stringsAsFactors = F)
chromosomes = c(1:24, "X", "Y")
gene_data = gene_data[gene_data$Chromosome.scaffold.name %in% chromosomes, ]
ngs = merge(ngs, unique(gene_data[,c("Gene.stable.ID", "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name", "Gene.type")]), 
            by.x = "Ensembl_ID_no_Dot", by.y = "Gene.stable.ID")

diff_genes = ngs[ngs$FDR <= FDR_cutoff & (ngs$FC <= 1 - fold_change_cutoff | ngs$FC >= 1 + fold_change_cutoff), ]
diff_genes$EKP_ID = NA

if(proteinOnly){
  diff_genes = diff_genes[diff_genes$Gene.type == "protein_coding", ]
}

# Use the Ensembl Gene Identifiers to map the genes to EKP identifiers
# Not all can be mapped. Here we count how many, and which
unmapped = c()
for(d in 1:nrow(diff_genes)){
  response = GET(paste0("http://lagavulin:9983/solr/collection1/select?q=term%3A%22", tolower(diff_genes$Ensembl_ID_no_Dot[d]), "%22+AND+semantictype%3A28+AND+knowledgebase%3A%22ensembl%22&wt=json"))
  response_content = fromJSON(content(response))
  identifiers = unlist(lapply(response_content$response$docs, function(x){return(x$gi)}))
  if(length(unique(identifiers)) == 1){
    diff_genes[d, "EKP_ID"] = unique(identifiers)
  } else {
    unmapped = append(unmapped, diff_genes$Ensembl_ID_no_Dot[d])
  }
}

diff_genes = diff_genes[!is.na(diff_genes$EKP_ID), ]
  
# Extract the triples from the knowledge graph
## Set up the database connection
triples = mongo(url = "mongodb://lagavulin:30017", db = "braindata", collection = "triples") 
concepts = mongo(url = "mongodb://lagavulin:30017", db = "braindata", collection = "concepts")
#EKP_genes = concepts$find(query = '{"st":"T028", "taxonomies":"homo sapiens"}')

rels = triples$find(query = paste0('{"s": {"$in" : ', toJSON(diff_genes$EKP_ID), '}, "o": {"$in" : ', toJSON(diff_genes$EKP_ID), '}}'))
rels = rels[rels$s != rels$o, ] # Remove selfref relationships
EKP_names = unique(data.frame(EKP_ID = c(rels$s, rels$o), EKP_name = c(rels$sn, rels$on)))

processes = triples$find(query = paste0('{"s": {"$in" : ', toJSON(diff_genes$EKP_ID), '}, "sco": "Physiology"}'))
processes2 = triples$find(query = paste0('{"o": {"$in" : ', toJSON(diff_genes$EKP_ID), '}, "scs": "Physiology"}'))

molecular_functions = concepts$find(query = '{"st": "T044"}', fields = '{"_id" : 1, "n" : 1}')
#physiological_functions = concepts$find(query = '{"st": "T039"}') # Niet helemaal zeker hierover
#neoplastic_process = concepts$find(query = '{"st": "T191"}') # Niet helemaal zeker hierover
#dysfunctions = concepts$find(query = '{"st": "T049"}') # Ook niet heel overtuigend
#cell_functions = concepts$find(query = '{"st": "T043"}') 
genetic_functions = concepts$find(query = '{"st": "T045"}', fields = '{"_id" : 1, "n" : 1}') 
#pathologic_functions = concepts$find(query = '{"st": "T046"}') 

processes = processes[processes$o %in% c(molecular_functions$'_id', genetic_functions$'_id'), ]
processes2 = processes2[processes2$s %in% c(molecular_functions$'_id', genetic_functions$'_id'), ]

all_processes = unique(data.frame(Gene = c(processes$s, processes2$o), Process = c(processes$o, processes2$s)))
process_frequency = table(all_processes$Process)
process_frequency = process_frequency[process_frequency > 1]

processes = processes[processes$o %in% names(process_frequency), ]
processes = processes[!is.na(processes$on), ]
processes2 = processes2[processes2$s %in% names(process_frequency), ]

# Create a process matrix -- ignore predicate information for now
process_graph = graph_from_edgelist(as.matrix(processes[,c("s", "on")]), directed = F)
process_matrix = as_adjacency_matrix(simplify(process_graph))
process_matrix = process_matrix[rownames(process_matrix) %in% as.character(processes$s), colnames(process_matrix) %in% processes$on]
rownames(process_matrix) = as.numeric(rownames(process_matrix))
colnames(process_matrix) = make.names(colnames(process_matrix))

# Shut down the connection
triples$disconnect()
concepts$disconnect()

#TODO: Nog alles simplifyen?
ppi = graph_from_edgelist(as.matrix(rels[,c("sn", "on")]))
ppi_undirected = graph_from_edgelist(as.matrix(rels[,c("sn", "on")]), directed = F)

Overall_Degree = degree(ppi, mode = "all")
eigen_values = eigen_centrality(ppi, directed = F)$vector
betweenness_values = betweenness(ppi, directed = F)
ppi_clusters = cluster_louvain(simplify(ppi_undirected))
ppi_clusters = data.frame(EKP_name = ppi_clusters$names, Cluster = ppi_clusters$membership)
ppi_clusters$Cluster = as.factor(ppi_clusters$Cluster)
#incoming = degree(ppi, mode = "in")
#outgoing = degree(ppi, mode = "out")

# Remove all genes that have only 1 edge
removal_genes = names(Overall_Degree[Overall_Degree == 1])
ppi = delete_vertices(ppi, removal_genes)

# Also remove these genes from the triples
rels = rels[which(!rels$sn %in% removal_genes & !rels$on %in% removal_genes), ]

# Specify the undirected predicates
undirected_predicates = c("interacts with", 
                          "forms protein complex with", 
                          "coexists with", 
                          "is compared with", 
                          "is functionally related to",
                          "does not interact with",
                          "does not coexist with",
                          "is associated with",
                          "binds with",
                          "is the same as",
                          "is spatially related to",
                          "ortholog is associated with")

outgoing = as.data.frame.matrix(table(rels$s[!rels$pn %in% undirected_predicates], rels$pn[!rels$pn %in% undirected_predicates]))
incoming = as.data.frame.matrix(table(rels$o[!rels$pn %in% undirected_predicates], rels$pn[!rels$pn %in% undirected_predicates]))
colnames(outgoing) = paste0(make.names(colnames(outgoing)), "_outgoing")
colnames(incoming) = paste0(make.names(colnames(incoming)), "_incoming")

undirected_rels = rels[rels$pn %in% undirected_predicates, c("s", "o", "pn")]
undirected_rels[,c("s", "o")] = t(apply(undirected_rels[,c("s", "o")], 1, sort))
undirected_rels = unique(undirected_rels)
undirected = merge(as.data.frame(table(undirected_rels$s, undirected_rels$pn)), as.data.frame(table(undirected_rels$o, undirected_rels$pn)), by = c("Var1", "Var2"), all = T)
undirected[is.na(undirected)] = 0
undirected$Total_undirected = undirected$Freq.x + undirected$Freq.y
undirected[,c("Freq.x", "Freq.y")] = NULL
undirected = spread(undirected, Var2, Total_undirected)
colnames(undirected) = make.names(colnames(undirected))

features = diff_genes[diff_genes$EKP_ID %in% c(rels$s, rels$o), c("Name", "FC", "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name", "EKP_ID")]
features = merge(features, EKP_names, by = "EKP_ID", all.x = T)
features = merge(features, incoming, by.x = "EKP_ID", by.y = "row.names", all.x = T)
features = merge(features, outgoing, by.x = "EKP_ID", by.y = "row.names", all.x = T)
features = merge(features, undirected, by.x = "EKP_ID", by.y = "Var1", all.x = T)
features = merge(features, as.data.frame(as.matrix(process_matrix)), by.x = "EKP_ID", by.y = "row.names", all.x = T)
features = merge(features, as.data.frame(Overall_Degree), by.x = "EKP_name", by.y = "row.names")
features = merge(features, as.data.frame(eigen_values), by.x = "EKP_name", by.y = "row.names")
features = merge(features, ppi_clusters, by = "EKP_name", all.x = T)
features[is.na(features)] = 0

refset = read.csv2("Reference sets/Genes from Farashi and Hazelett.csv", stringsAsFactors = F)
features$Class = factor(ifelse(features$Name %in% refset$x, "Target", "NonTarget"), levels = c("Target", "NonTarget"))
