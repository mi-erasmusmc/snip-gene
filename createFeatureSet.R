# Load the required packages
require(data.table)
require(rjson)

# Identify the negative genes
negatives = unique(unlist(lapply(combined$Candidate_genes, fromJSON)))
negatives = setdiff(negatives, combined$Target_Gene)

features = diff_genes[diff_genes$EKP_ID %in% c(rels$s, rels$o), c("Name", "Gene.start..bp.", "Gene.end..bp.", 
                                                                  "Chromosome.scaffold.name", "FC", "FDR", "diff_expression", "EKP_ID")]
features = merge(features, unique(data.frame(EKP_ID = c(rels$s, rels$o), EKP_Name = c(rels$sn, rels$on))), by = "EKP_ID", all.x = T)
features = merge(features, network_features, by.x = "EKP_Name", by.y = "Gene", all.x = T)
features = merge(features, pred_features, by = "EKP_ID", all.x = T)
features = merge(features, process_matrix, by.x = "EKP_ID", by.y = "row.names", all.x = T)

if(config$onlySelectNegatives){
  features = features[features$Name %in% c(as.character(combined$Target_Gene), negatives), ]
}
if(length(which(is.na(features))) > 0){
  features[is.na(features)] = 0
}
features$Class = factor(ifelse(features$Name %in% as.character(combined$Target_Gene), "Target", "NonTarget"), levels = c("Target", "NonTarget"))

fwrite(features, paste0("Raw data files/Complete feature set generated on ", Sys.Date(), ".csv"))
#refset = read.csv2("Reference sets/Genes from Farashi and Hazelett.csv", stringsAsFactors = F)