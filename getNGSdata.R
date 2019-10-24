# Load the required packages
require(httr)
require(rjson)

if(config$calculateAnew){
  # Select the relevant genes
  ngs = read.csv("GeneExpressionData/Galaxy37-[edgeR_DGE_on_2__design_matrix_prostate_unpaired.txt_-_differentially_expressed_genes].tabular.annotated.txt", 
                 stringsAsFactors = F, sep = "\t", row.names = 1)
  colnames(ngs) = c("Ensembl_ID", "Name", "log2FC", "logCPM", "LR", "p_value", "FDR")
  ngs$Ensembl_ID_no_Dot = sapply(ngs$Ensembl_ID, function(x){return(unlist(strsplit(x, "\\."))[1])})
  ngs$FC = 2^ngs$log2FC
  
  if(config$ExpressedGenesOnly){
    ngs = ngs[ngs$Ensembl_ID %in% expressed_genes, ]
  }
  
  # Get the locations of the genes
  gene_data = read.csv("GeneExpressionData/ENSEMBL gene locations downloaded from BioMart on 13-05-2019.txt", stringsAsFactors = F)
  chromosomes = c(1:22, "X", "Y")
  gene_data = gene_data[gene_data$Chromosome.scaffold.name %in% chromosomes, ]
  ngs = merge(ngs, unique(gene_data[,c("Gene.stable.ID", "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name", "Gene.type")]), 
              by.x = "Ensembl_ID_no_Dot", by.y = "Gene.stable.ID")
  
  ngs$EKP_ID = NA
  
  # Use the Ensembl Gene Identifiers to map the genes to EKP identifiers
  # Not all can be mapped. Here we count how many, and which
  unmapped = c()
  for(d in 1:nrow(ngs)){
    response = GET(paste0("http://lagavulin:9983/solr/collection1/select?q=term%3A%22", tolower(ngs$Ensembl_ID_no_Dot[d]), "%22+AND+semantictype%3A28+AND+knowledgebase%3A%22ensembl%22&wt=json"))
    response_content = fromJSON(content(response))
    identifiers = unlist(lapply(response_content$response$docs, function(x){return(x$gi)}))
    if(length(unique(identifiers)) == 1){
      ngs[d, "EKP_ID"] = unique(identifiers)
    } else {
      unmapped = append(unmapped, ngs$Ensembl_ID_no_Dot[d])
    }
  }
  ngs = ngs[!is.na(ngs$EKP_ID), ]
  
  ngs$expression = ifelse(ngs$FC > 1, "Upregulated", "Downregulated")
  ngs$expression[which(ngs$FC > 1 - config$fold_change_cutoff & ngs$FC < 1 + config$fold_change_cutoff)] = "Unchanged"
  
  write.csv2(ngs, paste0("Raw data files/Expressed genes and their differentiality ", todays_date, ".csv"), row.names = F)
} else {
  ngs = read.csv2(paste0("Raw data files/Expressed genes and their differentiality ", config$Data.date, ".csv"), stringsAsFactors = F)
}
  
if(config$diffGenesOnly){
  ngs = ngs[ngs$FDR <= config$FDR_cutoff & (ngs$FC <= 1 - config$fold_change_cutoff | ngs$FC >= 1 + config$fold_change_cutoff), ]
}

if(config$proteinOnly){
  ngs = ngs[ngs$Gene.type == "protein_coding", ]
}
