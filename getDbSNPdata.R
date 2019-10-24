# Get the SNP data from dbSNP
# Querying SNPs has to be done one by one due to an error in the rsnps package
# I therefore made this a one-time run script, instead of repeating it

# Load the required packages
require(xlsx)
require(DescTools)

# Read in the data from Farashi
farashi_xls = read.xlsx("Reference sets/41568_2018_87_MOESM1_ESM-3.xls", sheetName = "Supplementary TABLE 1", startRow = 2, endRow = 1288, stringsAsFactors = F)
farashi = farashi_xls[,c(1,3,4,5,7)]
farashi = farashi[!is.na(farashi$SNP.s.Genomic.Location) & !is.na(farashi$Target.assigned.e.Gene), ]

if(config$getSNPDataFromWeb){
  require(rsnps)
  
  # Remove entries withouth rs identifier
  farashi = farashi[-grep("chr", farashi$SNP.ID), ]
  
  snps = trimws(unique(farashi$SNP.ID))
  snps = gsub("_C", "", snps)
  snps = gsub("_A", "", snps)
  out = data.frame(Query = snps, Chromosome = NA, Marker = NA, Class = NA, Gene = NA, 
                   Alleles = NA, Major = NA, Minor = NA, MAF = NA, BP = NA, AncestralAllele = NA, 
                   stringsAsFactors = F)
  for(i in which(is.na(out$Marker))){
    print(i)
    snp_position = ncbi_snp_query(out$Query[i])
    if(nrow(snp_position) > 0){
      out[i,] = snp_position
    }
    Sys.sleep(2)
  }
  
  write.csv2(out, paste0("Raw data files/SNP info Farashi et al. Supplemental Materials 1 retrieved on ", todays_date, ".csv"), row.names = F)
} else {
  out = read.csv2(paste0("Raw data files//SNP info Farashi et al. Supplemental Materials 1 retrieved on ", config$Data.date, ".csv"), stringsAsFactors = F)
}

# Remove the SNPs that don't have a location
out = out[!is.na(out$BP), ]

# Modify farashi
farashi$Target.assigned.e.Gene = gsub("_(PSA)", "", farashi$Target.assigned.e.Gene)
farashi$Target.assigned.e.Gene = gsub(" (PSA)", "", farashi$Target.assigned.e.Gene, fixed = T)
farashi$Target.assigned.e.Gene = gsub(" (lncRNA)", "", farashi$Target.assigned.e.Gene, fixed = T)

# Update outdated gene identifiers
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene %in% c("MSMB1", "MSMB2")] = "MSMB"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "NCOA4-1"] = "NCOA4P1"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "NCOA4-3"] = "NCOA4P3"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "ANKRD5"] = "ANKEF1"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "C6orf228"] = "SMIM13"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "c-MYC, FoxA1 binding"] = "MYC"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "HoxB13"] = "HOXB13"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "LASS2"] = "CERS2"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "C10orf32"] = "BORCS7"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "LOC100505761"] = "RPARP-AS1"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "LOC100505495"] = "PCAT19"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "WDR52"] = "CFAP44"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "HCG4P6"] = "HCG4B"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "LOC285830"] = "HLA-F-AS1"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "RAB7L1"] = "RAB29"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "LOC284578"] = "MFSD4A-AS1"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "AGAP7"] = "AGAP7P"
farashi$Target.assigned.e.Gene[farashi$Target.assigned.e.Gene == "C2orf43"] = "LDAH"

split1 = strsplit(farashi$Target.assigned.e.Gene, ";")
farashi_unnested = farashi[rep(seq_len(nrow(farashi)), sapply(split1, length)), c("SNP.ID", "SNP.s.Genomic.Location", "Expeimental.approach", "reference")]

farashi_unnested$Target_Gene = unlist(split1)
farashi_unnested = unique(farashi_unnested)

split2 = strsplit(farashi_unnested$Target_Gene, ", ")
farashi_final = farashi_unnested[rep(seq_len(nrow(farashi_unnested)), sapply(split2, length)), c("SNP.ID", "SNP.s.Genomic.Location", "Expeimental.approach", "reference")]
farashi_final$Target_Gene = unlist(split2)
farashi_final = unique(farashi_final)

# Remove starting and trailing whitespaces from columns
farashi_final = as.data.frame(apply(farashi_final, 2, trimws))

# Remove entries which are based on a single eQTL study
farashi_final = farashi_final[!farashi_final$reference %in% c("Dadaev T. et al. 2018",
                                                              "Grisanzio, C. et al.  2012",
                                                              "X. xu et al. 2014",
                                                              "Thibodeau S.N.  et al. 2015"), ]

farashi_final = farashi_final[!farashi_final$SNP.s.Genomic.Location %in% c("Coding region",
                                                                           "exonic"), ]

farashi_final = merge(farashi_final, out, by.x = "SNP.ID", by.y = "Query")

# Combine SNP data with differentially expressed genes
for(i in 1:nrow(farashi_final)){
  chromosome_candidates = ngs[ngs$Chromosome.scaffold.name == farashi_final$Chromosome[i], ]
  chromosome_candidates$Overlap = Overlap(as.matrix(chromosome_candidates[,c("Gene.start..bp.", "Gene.end..bp.")]),
                                          c(farashi_final$BP[i] - config$delta_bp, farashi_final$BP[i] + config$delta_bp))
  farashi_final[i, "Candidate_genes"] = toJSON(unique(chromosome_candidates$Name[chromosome_candidates$Overlap > 0]))
  farashi_final[i, "nCandidates"] = length(unique(chromosome_candidates$Ensembl_ID[chromosome_candidates$Overlap > 0]))
  farashi_final[i, "TargetGeneInCandidates"] = ifelse(farashi_final$Target_Gene[i] %in% chromosome_candidates$Name[chromosome_candidates$Overlap > 0], T, F)
  farashi_final[i, "TargetGeneOnChromosome"] = ifelse(farashi_final$Target_Gene[i] %in% chromosome_candidates$Name, T, F)
}

write.csv2(farashi_final, paste0("Raw data files/Genes and gene candidates identified on ", todays_date, ".csv"), row.names = F)
