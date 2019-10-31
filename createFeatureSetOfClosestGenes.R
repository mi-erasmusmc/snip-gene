# Create a training set that consists of the closest genes as positives and negatives, and a test set of the other gene candidates

# Load the required packages
require(data.table)
require(rjson)

farashi_final = farashi_final[!is.na(farashi_final$TargetGeneOnChromosome), ]

# Identify the positive genes for the training set
positives_train = farashi_final[farashi_final$DistanceRankTargetGene == 1, ]

# Identify the negative genes for the training set
negatives = farashi_final[farashi_final$DistanceRankTargetGene > 1, ]
negatives_train = lapply(na.omit(negatives$CandidatesSortedOnDistance), fromJSON)
negatives_train = sapply(negatives_train, function(x){x[1]})

features = ngs[ngs$EKP_ID %in% c(rels$s, rels$o), c("Name", "Gene.start..bp.", "Gene.end..bp.", 
                                                    "Chromosome.scaffold.name", "FC", "FDR", "expression", "EKP_ID")]

features = merge(features, unique(data.frame(EKP_ID = c(rels$s, rels$o), EKP_Name = c(rels$sn, rels$on))), by = "EKP_ID", all.x = T)
features = merge(features, pred_features, by = "EKP_ID", all.x = T)
features = merge(features, process_matrix, by.x = "EKP_ID", by.y = "row.names", all.x = T)

features_train = features[features$Name %in% c(as.character(positives_train$Target_Gene), negatives_train), ]

if(length(which(is.na(features_train))) > 0){
  features_train[is.na(features_train)] = 0
}
features_train$Class = factor(ifelse(features_train$Name %in% as.character(positives_train$Target_Gene), "Target", "NonTarget"), levels = c("Target", "NonTarget"))

fwrite(features_train, paste0("Raw data files/First-gene training set generated on ", todays_date, ".csv"))
#refset = read.csv2("Reference sets/Genes from Farashi and Hazelett.csv", stringsAsFactors = F)

# Create the test-set
negatives_test = lapply(na.omit(negatives$CandidatesSortedOnDistance), fromJSON)
negatives_test = unlist(sapply(negatives_test, function(x){x[2:length(x)]}))

features = features[features$Name %in% as.character(negatives_test), ]

if(length(which(is.na(features))) > 0){
  features[is.na(features)] = 0
}

features$Class = factor(ifelse(features$Name %in% as.character(negatives$Target_Gene), "Target", "NonTarget"), levels = c("Target", "NonTarget"))

fwrite(features, paste0("Raw data files/First-gene test set generated on ", todays_date, ".csv"))