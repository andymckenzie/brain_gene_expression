library(edgeR)
library(DGCA)
library(limma)
library(HGNChelper)


setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

############################
# functions

#https://www.biostars.org/p/72846/
cpm_tmm <- function(counts){
    d = DGEList(counts = counts)
    d = calcNormFactors(d, method = "TMM")
    return(cpm(d, normalized.lib.sizes = TRUE))
}

###################################
# tasic

metadata = read.csv("data/tasic_cell_metadata.csv")
counts = read.csv("data/tasic_genes_counts.csv")
rownames(counts) = counts[ , 1]
counts = counts[ , -1]
tasic_counts = counts
tasic_norm = cpm_tmm(counts)

sub_class = table(metadata$sub_class)
sub_class = sub_class[sub_class >= 15]
sub_class = sub_class[!names(sub_class) == ""]
cell_types = names(sub_class)

tasic_cell_types = metadata$sub_class
tasic_cell_types[grepl("L", tasic_cell_types)] = "Neuron"
tasic_cell_types[grepl("Sst", tasic_cell_types)] = "Neuron"
tasic_cell_types = gsub("Vip", "Neuron", tasic_cell_types)
tasic_cell_types = gsub("Th", "Neuron", tasic_cell_types)
tasic_cell_types = gsub("Pvalb", "Neuron", tasic_cell_types)
tasic_cell_types = gsub("Ndnf", "Neuron", tasic_cell_types)
tasic_cell_types = make.names(tasic_cell_types)

tasic_norm_log = log(tasic_norm + 1, 2)

##########################
#convert to cleaned human gene sybmols where possible
source("/Users/amckenz/Documents/github/alzolig/convert_mgi_to_hgnc.R")
gene_symbols = convert_mgi_to_hgnc(rownames(tasic_counts), genome = "mm10", return_df = FALSE)
gene_symbols = switchGenesToHGCN(gene_symbols)
rownames(tasic_counts) = make.unique(gene_symbols)
collapsed_df = collapseRows(tasic_counts, rowGroup = gene_symbols, rowID = make.unique(gene_symbols), method="MaxMean")
tasic_counts = tasic_counts[collapsed_df$selectedRow, ]
tasic_norm_log = tasic_norm_log[collapsed_df$selectedRow, ]
gene_symbols = gene_symbols[collapsed_df$selectedRow]
rownames(tasic_counts) = gene_symbols
#remove the spike-ins
mt_spike_in = grepl("MT_", rownames(tasic_counts), fixed = TRUE)
tasic_counts = tasic_counts[!mt_spike_in, ]
tasic_norm_log = tasic_norm_log[!mt_spike_in, ]

#############################
# dfxp

dfxp_function <- function(cell_type, list_other_cells){

  expr_average = rowMeans(as.matrix(tasic_norm_log[ , tasic_cell_types == cell_type]))

  celltypes_contrast = tasic_cell_types
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)

  tasic_voom_res = voom(tasic_counts, design, plot = TRUE)
	fit = lmFit(tasic_voom_res, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
	fit2 = eBayes(fit2)
	toptable = topTable(fit2, adjust = "BH", sort.by = "none", number = nrow(fit2), coef = 1, confint = TRUE)

  toptable$expr_average = expr_average
  source("shrink_logFC.R")
  toptable$fc_zscore = shrink_pval_and_avg_expr(toptable)
  toptable$genes = rownames(toptable)
	toptable = toptable[order(toptable$fc_zscore, decreasing = TRUE), ]

	return(toptable)

}

tasic_oli_df = dfxp_function("Oligodendrocyte", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))
tasic_neu_df = dfxp_function("Neuron", list("Oligodendrocyte", "Astrocyte", "Microglia", "Endothelial", "OPC"))
tasic_ast_df = dfxp_function("Astrocyte", list("Neuron", "Oligodendrocyte", "Microglia", "Endothelial", "OPC"))
tasic_mic_df = dfxp_function("Microglia", list("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial", "OPC"))
tasic_end_df = dfxp_function("Endothelial", list("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC"))
tasic_opc_df = dfxp_function("OPC", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))

saveRDS(tasic_oli_df, "data/dfxp/tasic_oli_df_dfxp.rds")
saveRDS(tasic_neu_df, "data/dfxp/tasic_neu_df_dfxp.rds")
saveRDS(tasic_ast_df, "data/dfxp/tasic_ast_df_dfxp.rds")
saveRDS(tasic_mic_df, "data/dfxp/tasic_mic_df_dfxp.rds")
saveRDS(tasic_end_df, "data/dfxp/tasic_end_df_dfxp.rds")
saveRDS(tasic_opc_df, "data/dfxp/tasic_opc_df_dfxp.rds")
