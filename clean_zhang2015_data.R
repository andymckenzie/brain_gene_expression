library(XLConnect)
library(limma)
library(DGCA)
library(HGNChelper)


setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

############################
# zhang 2015 mouse data

#downloaded the "supplementary excel files" for each sample from
#http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52564
filenames = list.files("data/GSE52564_RAW", pattern="*.xls", full.names=TRUE)
#may get Error: OutOfMemoryError (Java): Java heap space due to XLConnect here
ldf = lapply(filenames, readWorksheetFromFile, sheet = 1)
barres_merged = Reduce(function(...) merge(..., by = "gene.symbol", all = T),
	ldf)
row.names(barres_merged) = barres_merged[ , 1]
barres_merged = barres_merged[ , -1]

#following Gordon Smyth's advice for best practices
#https://support.bioconductor.org/p/56275/

#log2 transform the data
#fpkm values are already normalized
barres_merged_log = log(barres_merged + 1, 2)
#order from GEO GSM1269903 Astrocyte1 GSM1269904 Astrocyte2 GSM1269905 Neuron1 GSM1269906 Neuron2 GSM1269907 OPC1 GSM1269908 OPC2 GSM1269909 NFO1 GSM1269910 NFO2 GSM1269911 MO1 GSM1269912 MO2 GSM1269913 Microglia1 GSM1269914 Microglia2 GSM1269915 Endothelial1 GSM1269916 Endothelial2 GSM1269917 WC1 GSM1269918 WC2 GSM1269919 WC3
sample_types = c("Ast", "Ast", "Neu", "Neu", "Opc", "Opc",
	"Nfo", "Nfo", "Mol", "Mol", "Mic", "Mic", "Endo", "Endo",
	"WC", "WC", "WC")

##########################
#convert to cleaned human gene sybmols where possible
source("/Users/amckenz/Documents/github/alzolig/convert_mgi_to_hgnc.R")
gene_symbols = convert_mgi_to_hgnc(rownames(barres_merged_log), genome = "mm10", return_df = FALSE)
gene_symbols = switchGenesToHGCN(gene_symbols)
rownames(barres_merged_log) = make.unique(gene_symbols)
collapsed_df = collapseRows(barres_merged_log, rowGroup = gene_symbols, rowID = make.unique(gene_symbols), method="MaxMean")
barres_merged_log = barres_merged_log[collapsed_df$selectedRow, ]
gene_symbols = gene_symbols[collapsed_df$selectedRow]
rownames(barres_merged_log) = gene_symbols

####################
#dfxp

dfxp_function <- function(cell_type, list_other_cells){

	expr_average = rowMeans(barres_merged_log[ , which(sample_types %in% cell_type)])

  celltypes_contrast = sample_types
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
	fit = lmFit(barres_merged_log, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
	fit2 = eBayes(fit2, trend = TRUE)
	toptable = topTable(fit2, adjust = "BH", sort.by = "none", number = nrow(fit2), coef = 1, confint = TRUE)

  toptable$expr_average = expr_average
  source("shrink_logFC.R")
  toptable$fc_zscore = shrink_pval_and_avg_expr(toptable)
  toptable$genes = rownames(toptable)
	toptable = toptable[order(toptable$fc_zscore, decreasing = TRUE), ]

	return(toptable)

}

z15_oli_df = dfxp_function("Mol", list("Neu", "Ast", "Mic", "Endo"))
z15_neu_df = dfxp_function("Neu", list("Mol", "Ast", "Mic", "Endo", "Opc"))
z15_ast_df = dfxp_function("Ast", list("Neu", "Mol", "Mic", "Endo", "Opc"))
z15_mic_df = dfxp_function("Mic", list("Neu", "Ast", "Mol", "Endo", "Opc"))
z15_end_df = dfxp_function("Endo", list("Neu", "Ast", "Mic", "Mol", "Opc"))
z15_opc_df = dfxp_function("Opc", list("Neu", "Ast", "Mic", "Endo"))

saveRDS(z15_oli_df, "data/dfxp/z15_oli_df_dfxp.rds")
saveRDS(z15_neu_df, "data/dfxp/z15_neu_df_dfxp.rds")
saveRDS(z15_ast_df, "data/dfxp/z15_ast_df_dfxp.rds")
saveRDS(z15_mic_df, "data/dfxp/z15_mic_df_dfxp.rds")
saveRDS(z15_end_df, "data/dfxp/z15_end_df_dfxp.rds")
saveRDS(z15_opc_df, "data/dfxp/z15_opc_df_dfxp.rds")
