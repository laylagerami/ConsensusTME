# Install package
devtools::install_github("cansysbio/ConsensusTME")

# Load libary
library(ConsensusTME)

# Read expression matrix
bulkExpMatrix <- as.matrix(read.delim(bulkGeneExpression.txt, row.names = 1))

# Analysis for Ovarian Cancer with ssGSEA
ConsensusTME::consensusTMEAnalysis(bulkExpMatrix, cancer = "OV", statMethod = "ssgsea")