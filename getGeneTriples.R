# Get the gene triples

# Load the required packages
require(mongolite)
require(igraph)
require(data.table)
require(rjson)

if(config$getFromKG){
  # Extract the triples from the knowledge graph
  ## Set up the database connection
  triples = mongo(url = "mongodb://lagavulin:30017", db = "braindata", collection = "triples") 
  concepts = mongo(url = "mongodb://lagavulin:30017", db = "braindata", collection = "concepts")
  #EKP_genes = concepts$find(query = '{"st":"T028", "taxonomies":"homo sapiens"}')
  
  rels = triples$find(query = paste0('{"s": {"$in" : ', toJSON(diff_genes$EKP_ID), '}, "o": {"$in" : ', toJSON(diff_genes$EKP_ID), '}}'))
  rels = merge(rels, diff_genes[,c("EKP_ID", "diff_expression")], by.x = "o", by.y = "EKP_ID")
  rels = merge(rels, diff_genes[,c("EKP_ID", "diff_expression")], by.x = "s", by.y = "EKP_ID", suffixes = c("_Object", "_Subject"))
  rels$incoming = paste0(rels$diff_expression_Subject, "_", rels$pn)
  rels$outgoing = paste0(rels$pn, "_", rels$diff_expression_Object)
  
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
  
  # Only include processes that are associated with multiple genes
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
  if(config$saveTriples){
    rels$a = NULL #Eruit gehaald vanwege de write, maar misschien in de toekomst erin laten?
    rels$m = NULL #Eruit gehaald vanwege de write, maar misschien in de toekomst erin laten?
    
    fwrite(rels, file = paste0("Raw data files/Triples extracted from the knowledge graph on ", todays_date, ".csv"), sep = ";")
    fwrite(as.data.frame(as.matrix(process_matrix)), paste0("Raw data files/Process matrix extracted from the knowledge graph on ", todays_date, ".csv"), sep = ";")
  }
} else {
  rels = as.data.frame(fread(paste0("Raw data files/Triples extracted from the knowledge graph on ", config$Data.date, ".csv")))
  process_matrix = as.data.frame(fread(paste0("Raw data files/Process matrix extracted from the knowledge graph on", config$Data.date, ".csv")))
}