# Calculate which proteins are expressed within the normal and the prostate cancer data
print("Identifying expressed genes")

Expression_percentage = 0.5
nReads_cutoff = 2

reads = read.delim("GeneExpressionData/expression_matrix_prostate_clean.txt", stringsAsFactors = F, row.names = 1)

cols = read.delim("GeneExpressionData/design_matrix_prostate_unpaired.txt", stringsAsFactors = F)
cols$name = make.names(cols$samplename)

# Split out the table into two pieces
normal = reads[, colnames(reads) %in% cols$name[cols$condition == "normal"]]
cancer = reads[, colnames(reads) %in% cols$name[cols$condition == "cancer"]]

# Create binary matrices to identify in how many patients reads are found
normal_binary = normal
normal_binary[normal_binary < nReads_cutoff] = F
normal_binary[normal_binary >= nReads_cutoff] = T
normal_binary$total = rowSums(normal_binary)
normal_binary = normal_binary[normal_binary$total > (ncol(normal_binary) - 1) * Expression_percentage,]

cancer_binary = cancer
cancer_binary[cancer_binary < nReads_cutoff] = F
cancer_binary[cancer_binary >= nReads_cutoff] = T
cancer_binary$total = rowSums(cancer_binary)
cancer_binary = cancer_binary[cancer_binary$total > (ncol(cancer_binary) - 1) * Expression_percentage,]

expressed_genes = union(row.names(cancer_binary), row.names(normal_binary))
