library(XLConnect)
library(limma)
library(DGCA)

setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

cutoff_thresh = 0.5

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
#log2(x + 0.25) is added by edgeR for stabilization https://www.biostars.org/p/102568/
#fpkm values are already normalized
barres_merged_log = log(barres_merged + 1, 2)
#order from GEO GSM1269903 Astrocyte1 GSM1269904 Astrocyte2 GSM1269905 Neuron1 GSM1269906 Neuron2 GSM1269907 OPC1 GSM1269908 OPC2 GSM1269909 NFO1 GSM1269910 NFO2 GSM1269911 MO1 GSM1269912 MO2 GSM1269913 Microglia1 GSM1269914 Microglia2 GSM1269915 Endothelial1 GSM1269916 Endothelial2 GSM1269917 WC1 GSM1269918 WC2 GSM1269919 WC3
sample_types = c("Ast", "Ast", "Neu", "Neu", "Opc", "Opc",
	"Nfo", "Nfo", "Mol", "Mol", "Mic", "Mic", "Endo", "Endo",
	"WC", "WC", "WC")

str(barres_merged_log)

dfxp_function <- function(cell_type, list_other_cells){

  #for a given cell type that you are considering
	#remove that gene if it has too low of expression *in that cell type*
	#find the average of each gene in its own cell type
	expr_average = rowMeans(barres_merged_log[ , which(sample_types %in%
		cell_type)])
  str(expr_average)
  print(sum(is.na(expr_average)))
	cutoff = as.numeric(quantile(expr_average, cutoff_thresh, names = FALSE))
	barres_trimmed = barres_merged_log[(expr_average > cutoff), ]

  celltypes_contrast = sample_types
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
	fit = lmFit(barres_trimmed, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
	fit2 = eBayes(fit2, 0.01, trend = TRUE)
  print(head(fit2$coefficients))
	tT = topTable(fit2, adjust = "BH", sort.by="B", number = nrow(fit2), coef = 1, confint = TRUE)

  #calculate the average vs all of the other cell types
	for(i in 1:length(list_other_cells)){
    mean_other = rowMeans(barres_trimmed[ , which(sample_types %in% list_other_cells[[i]])])
		if(i == 1){
			mean_others = data.frame(mean_other)
		} else {
			mean_others = cbind(mean_others, mean_other)
		}
    str(mean_others)
	}
  mean_fc = rowMeans(barres_trimmed[ , which(sample_types %in% cell_type)]) /
      rowMeans(mean_others)
	mean_fc_names = data.frame(row.names(barres_trimmed), mean_fc)
	names(mean_fc_names) = c("barres_names", "mean_fc")

	#since already sorted, need to merge via rownames
	toptable = merge(tT, mean_fc_names, by.x = "row.names", by.y = "barres_names")
	toptable = toptable[order(toptable$t, decreasing = TRUE), ]

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
