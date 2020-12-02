# Install package
devtools::install_github("cansysbio/ConsensusTME")

# Load libary
library(ConsensusTME)

# Read expression matrix
bulkExpMatrix <- as.matrix(read.delim("../data-raw/CIBERSORT_LM22.txt", row.names = 1))

# Analysis for Ovarian Cancer with ssGSEA
res = ConsensusTME::consensusTMEAnalysis(bulkExpMatrix, cancer = "OV", statMethod = "ssgsea")

# Print results
head(res)

# Create heatmap
library(ComplexHeatmap)

# Function to save heatmap image
saveHeatmap<-function(data){ 
  png( file=enc2native("figs/OC_heatmap.png"),width=800, height=650) 
  print(Heatmap(data,
                name = "NES (Normalised Enrichment Score)", #title of legend
                column_title = "Cell Types 2", row_title = "Cell Types 1",
                row_names_gp = gpar(fontsize = 7)
                ) 
  )
  dev.off() 
} 

saveHeatmap(res)

### REPEAT ANALYSIS FOR DIFFERENT CANCER TYPE
# Analysis for Ovarian Cancer with ssGSEA
res = ConsensusTME::consensusTMEAnalysis(bulkExpMatrix, cancer = "BRCA", statMethod = "ssgsea")

# Print results
head(res)

# Function to save heatmap image
saveHeatmap<-function(data){ 
  png( file=enc2native("figs/BRCA_heatmap.png"),width=800, height=650) 
  print(Heatmap(data,
                name = "NES (Normalised Enrichment Score)", #title of legend
                column_title = "Cell Types 2", row_title = "Cell Types 1",
                row_names_gp = gpar(fontsize = 7)
  ) 
  )
  dev.off() 
} 

saveHeatmap(res)

### PATHWAY ENRICHMENT
head(bulkExpMatrix)

# Extract Gene expression for Naive B cells
NBC_data = data.frame(bulkExpMatrix[,1])
colnames(NBC_data) = "Naive B Cell Expression"

# First we need to convert HGNC to Entrez ID
library('org.Hs.eg.db')
NBC_data$Entrez = mapIds(org.Hs.eg.db, 
       rownames(NBC_data),
       'ENTREZID', 
       'SYMBOL')

# Now we can carry out enrichment :)
library(clusterProfiler)

geneList = as.numeric(as.character(NBC_data$`Naive B Cell Expression`))
names(geneList) = as.character(NBC_data$Entrez)
geneList = geneList[!is.na(names(geneList))]
geneList = sort(geneList, decreasing = T)
nbc_gsea = gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 nPerm        = 1000,
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.5,
                 verbose      = T)

