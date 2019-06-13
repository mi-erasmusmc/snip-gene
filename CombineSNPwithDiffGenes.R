# Combine SNP data with differentially expressed genes
require(DescTools)

combined = farashi_final

if(config$diffGenesOnly){
  for(i in 1:nrow(combined)){
    chromosome_candidates = diff_genes[diff_genes$Chromosome.scaffold.name == combined$Chromosome[i], ]
    chromosome_candidates$Overlap = Overlap(as.matrix(chromosome_candidates[,c("Gene.start..bp.", "Gene.end..bp.")]),
                                           c(combined$BP[i] - config$delta_bp, combined$BP[i] + config$delta_bp))
    combined[i, "Candidate_genes"] = toJSON(unique(chromosome_candidates$Name[chromosome_candidates$Overlap > 0]))
    combined[i, "nCandidates"] = length(unique(chromosome_candidates$Ensembl_ID[chromosome_candidates$Overlap > 0]))
    combined[i, "TargetGeneInCandidates"] = ifelse(combined$Target_Gene[i] %in% chromosome_candidates$Name[chromosome_candidates$Overlap > 0], T, F)
    combined[i, "TargetGeneOnChromosome"] = ifelse(combined$Target_Gene[i] %in% chromosome_candidates$Name, T, F)
    }
} else {
  for(i in 1:nrow(combined)){
    chromosome_candidates = ngs[ngs$Chromosome.scaffold.name == combined$Chromosome[i], ]
    chromosome_candidates$Overlap = Overlap(as.matrix(chromosome_candidates[,c("Gene.start..bp.", "Gene.end..bp.")]),
                                            c(combined$BP[i] - config$delta_bp, combined$BP[i] + config$delta_bp))
    combined[i, "Candidate_genes"] = toJSON(unique(chromosome_candidates$Name[chromosome_candidates$Overlap > 0]))
    combined[i, "nCandidates"] = length(unique(chromosome_candidates$Ensembl_ID[chromosome_candidates$Overlap > 0]))
    combined[i, "TargetGeneInCandidates"] = ifelse(combined$Target_Gene[i] %in% chromosome_candidates$Name[chromosome_candidates$Overlap > 0], T, F)
    combined[i, "TargetGeneOnChromosome"] = ifelse(combined$Target_Gene[i] %in% chromosome_candidates$Name, T, F)
  }
}