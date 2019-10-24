# Create the features that are based on the direct and the predicates from the direct and the indirect relationships
require(data.table)

if(config$generatePredicateFeatures){
  print("Creating the predicate feature set")
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
  
  # Direct relationships
  print("Creating features for direct relationships")
  outgoing = as.data.frame.matrix(table(rels$s[!rels$pn %in% undirected_predicates], rels$outgoing[!rels$pn %in% undirected_predicates]))
  incoming = as.data.frame.matrix(table(rels$o[!rels$pn %in% undirected_predicates], rels$incoming[!rels$pn %in% undirected_predicates]))

  undirected_rels = rels[rels$pn %in% undirected_predicates, c("s", "o")]
  undirected_rels = as.data.frame(t(apply(undirected_rels, 1, sort)))
  undirected_rels$pn = rels[rels$pn %in% undirected_predicates, "pn"]
  undirected_rels = unique(undirected_rels)
  
  # Add the differential expression of the proteins
  undirected_rels = merge(undirected_rels, ngs[,c("EKP_ID", "expression")], by.x = "V1", by.y = "EKP_ID")
  undirected_rels = merge(undirected_rels, ngs[,c("EKP_ID", "expression")], by.x = "V2", by.y = "EKP_ID", suffixes = c("_V1", "_V2"))
  
  undirected_with_expression = data.frame(Gene = c(undirected_rels$V1, undirected_rels$V2), 
                                          predicate = c(paste0(undirected_rels$pn, "_", undirected_rels$expression_V2),
                                                        paste0(undirected_rels$pn, "_", undirected_rels$expression_V1)))
  
  undirected = as.data.frame.matrix(table(undirected_with_expression$Gene, undirected_with_expression$predicate))

  direct = merge(incoming, outgoing, by = "row.names", all = T)
  direct = merge(undirected, direct, by.x = "row.names", by.y = "Row.names", all = T)
  colnames(direct) = make.names(colnames(direct))
  
  # Clean up
  incoming = NULL
  outgoing = NULL
  undirected = NULL
  undirected_with_expression = NULL
  
  # Indirect relationships
  ## Directed predicates only
  print("Creating features for indirect relationships")
  indirect_dd = merge(as.data.table(rels[-which(rels$pn %in% undirected_predicates), c("s", "pn", "o", "expression_Subject", "expression_Object")]), 
                   as.data.table(rels[-which(rels$pn %in% undirected_predicates), c("s", "pn", "o", "expression_Object")]), 
                   by.x = "o", by.y = "s", allow.cartesian = T)
  indirect_dd = unique(indirect_dd)
  colnames(indirect_dd) = c("Intermediate", "Subject", "Pred1", "expression_Subject", "expression_Intermediate", "Pred2", "Object", "expression_Object")
  
  ## Also include the undirected predicates
  indirect_du1 = merge(as.data.table(rels[-which(rels$pn %in% undirected_predicates), c("s", "pn", "o", "expression_Subject", "expression_Object")]), 
                      as.data.table(rels[rels$pn %in% undirected_predicates, c("s", "pn", "o", "expression_Object")]), 
                      by.x = "o", by.y = "s", allow.cartesian = T)
  colnames(indirect_du1) = c("Intermediate", "Subject", "Pred1", "expression_Subject", "expression_Intermediate", "Pred2", "Object", "expression_Object")
  
  indirect_du2 = merge(as.data.table(rels[-which(rels$pn %in% undirected_predicates), c("s", "pn", "o", "expression_Subject", "expression_Object")]), 
                       as.data.table(rels[rels$pn %in% undirected_predicates, c("s", "pn", "o", "expression_Subject")]), 
                       by = "o", allow.cartesian = T)
  colnames(indirect_du2) = c("Intermediate", "Subject", "Pred1", "expression_Subject", "expression_Intermediate", "Object", "Pred2", "expression_Object")
  
  indirect_du = unique(rbind(indirect_du1, indirect_du2))
  
  
  # Also other way around
  indirect_ud1 = merge(as.data.table(rels[rels$pn %in% undirected_predicates, c("s", "pn", "o", "expression_Subject", "expression_Object")]),
                       as.data.table(rels[-which(rels$pn %in% undirected_predicates), c("s", "pn", "o", "expression_Object")]), 
                       by.x = "o", by.y = "s", allow.cartesian = T)
  colnames(indirect_ud1) = c("Intermediate", "Subject", "Pred1", "expression_Subject", "expression_Intermediate", "Pred2", "Object", "expression_Object")
  
  indirect_ud2 = merge(as.data.table(rels[rels$pn %in% undirected_predicates, c("s", "pn", "o", "expression_Object")]),
                       as.data.table(rels[-which(rels$pn %in% undirected_predicates), c("s", "pn", "o", "expression_Subject", "expression_Object")]), 
                       by = "s", allow.cartesian = T)
  colnames(indirect_ud2) = c("Intermediate", "Pred1", "Subject", "expression_Subject", "Pred2", "Object", "expression_Intermediate", "expression_Object")
  
  indirect_ud = unique(rbind(indirect_ud1, indirect_ud2))
  
  # Combineer alle resultaten
  indirect = rbind(indirect_dd, indirect_du, indirect_ud)
  indirect = unique(indirect[indirect$Subject != indirect$Object, ])
  
  indirect$outgoing = paste0(indirect$Pred1, "_", indirect$expression_Intermediate, "_", indirect$Pred2, "_", indirect$expression_Object)
  indirect$incoming = paste0(indirect$expression_Subject, "_", indirect$Pred1, "_", indirect$expression_Intermediate, "_", indirect$Pred2)
  
  # Clean up
  indirect_dd = NULL
  indirect_du = NULL
  indirect_ud = NULL
  
  # Finalize
  indirect_outgoing= as.data.frame.matrix(table(indirect$Subject, indirect$outgoing))
  indirect_incoming = as.data.frame.matrix(table(indirect$Object, indirect$incoming))

  indirect = merge(indirect_incoming, indirect_outgoing, by = "row.names", all = T)
  colnames(indirect) = make.names(colnames(indirect))
  
  pred_features = merge(direct, indirect, by = "Row.names", all = T)
  colnames(pred_features)[1] = "EKP_ID"
  
  pred_features[is.na(pred_features)] = 0
  
  fwrite(pred_features, paste0("Raw data files/Predicate features generated on ", todays_date, ".csv"))
} else {
  print("Reading in the previously created predicate features")
  pred_features = as.data.frame(fread(paste0("Raw data files/Predicate features generated on ", config$Data.date, ".csv")))
}