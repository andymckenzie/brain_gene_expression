library(readxl)
library(edgeR)
library(DGCA)
library(limma)
library(HGNChelper)

setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

################################
# sharma
#read in the sharma RPKM data set
sharma = read_excel("data/sharma_rna_nn.4160-S5.xlsx")
sharma_num = sharma[ , !colnames(sharma) %in% "GeneName"]
#NAs are 0s in this data set
sharma_num[is.na(sharma_num)] = 0
# > sum(is.na(sharma_num))
# [1] 0
# > sum(sharma_num == 0, na.rm = TRUE)
# [1] 97127
#this plot shows that normalizing the sharma data is not necessary..
#plot(colMeans(as.matrix(sharma_num)), colMedians(as.matrix(sharma_num)))
rownames(sharma_num) = make.unique(sharma[ , colnames(sharma) %in% "GeneName"])

sharma_cell_types = colnames(sharma_num)
sharma_cell_types[grepl("Oligodendrocytes div1", sharma_cell_types)] = "OPC"
sharma_cell_types[grepl("Oligodendrocytes div4", sharma_cell_types)] = "Oligodendrocyte"
sharma_cell_types[grepl("adult microglia", sharma_cell_types)] = "Microglia"
sharma_cell_types[grepl("cortical neurons div10", sharma_cell_types)] = "Neuron"
sharma_cell_types[grepl("Astrocytes", sharma_cell_types)] = "Astrocyte"
sharma_cell_types = make.names(sharma_cell_types)

sharma_num_log = log(sharma_num + 1, 2)

##########################
#convert to cleaned human gene sybmols where possible
source("/Users/amckenz/Documents/github/alzolig/convert_mgi_to_hgnc.R")
#sum(toupper(sharma[ , colnames(sharma) %in% "GeneName"]) != gene_symbols) ; 1699
gene_symbols = convert_mgi_to_hgnc(sharma[ , colnames(sharma) %in% "GeneName"], genome = "mm10", return_df = FALSE)
gene_symbols = switchGenesToHGCN(gene_symbols)
rownames(sharma_num_log) = make.unique(gene_symbols)
collapsed_df = collapseRows(sharma_num_log, rowGroup = gene_symbols, rowID = make.unique(gene_symbols), method = "MaxMean")
sharma_num_log = sharma_num_log[collapsed_df$selectedRow, ]
gene_symbols = gene_symbols[collapsed_df$selectedRow]
rownames(sharma_num_log) = gene_symbols

####################
#dfxp

dfxp_function <- function(cell_type, list_other_cells){

  expr_average = rowMeans(as.matrix(sharma_num_log[ , sharma_cell_types == cell_type]))

  celltypes_contrast = sharma_cell_types
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
	fit = lmFit(sharma_num_log, design)
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

sharma_oli_df = dfxp_function("Oligodendrocyte", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))
sharma_neu_df = dfxp_function("Neuron", list("Oligodendrocyte", "Astrocyte", "Microglia", "Endothelial", "OPC"))
sharma_ast_df = dfxp_function("Astrocyte", list("Neuron", "Oligodendrocyte", "Microglia", "Endothelial", "OPC"))
sharma_mic_df = dfxp_function("Microglia", list("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial", "OPC"))
sharma_opc_df = dfxp_function("OPC", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))

saveRDS(sharma_oli_df, "data/dfxp/sharma_oli_df_dfxp.rds")
saveRDS(sharma_neu_df, "data/dfxp/sharma_neu_df_dfxp.rds")
saveRDS(sharma_ast_df, "data/dfxp/sharma_ast_df_dfxp.rds")
saveRDS(sharma_mic_df, "data/dfxp/sharma_mic_df_dfxp.rds")
saveRDS(sharma_opc_df, "data/dfxp/sharma_opc_df_dfxp.rds")
