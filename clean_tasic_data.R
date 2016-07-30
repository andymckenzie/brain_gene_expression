library(edgeR)
library(DGCA)

cutoff_thresh_percentile = 0.5

setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

############################
# functions

std_error <- function(x) sd(x)/sqrt(length(x))

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

# colnames(tasic_total_df)[grepl("Sst", colnames(tasic_total_df))] =
#   c("Sst Cbln4_mean", "Sst Cbln4_se", "Sst Cbln4_mean_log", "Sst Cbln4_se_log",
#   "Sst Chodl_mean", "Sst Chodl_se", "Sst Chodl_mean_log", "Sst Chodl_se_log",
#   "Sst Nr2f2_mean", "Sst Nr2f2_se", "Sst Nr2f2_mean_log", "Sst Nr2f2_se_log")

tasic_norm_log = log(tasic_norm + 1, 2)

dfxp_function <- function(cell_type, list_other_cells){

	#for a given cell type that you are considering
	#remove that gene if it has too low of expression *in that cell type*
	#find the average of each gene in its own cell type
  expr_average = rowMeans(as.matrix(tasic_norm_log[ , tasic_cell_types == cell_type]))
	cutoff = as.numeric(quantile(expr_average, cutoff_thresh_percentile, names = FALSE))
	tasic_trimmed = tasic_norm_log[(expr_average > cutoff), ]

  celltypes_contrast = tasic_cell_types
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
	fit = lmFit(tasic_trimmed, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
	fit2 = eBayes(fit2, 0.01, trend = TRUE)
  print(head(fit2$coefficients))
	tT = topTable(fit2, adjust = "BH", sort.by="B", number = nrow(fit2), coef = 1, confint = TRUE)
	# tT = tT[tT$P.Value < pval_thresh, ]

  #calculate the average vs all of the other cell types
	for(i in 1:length(list_other_cells)){
    mean_other = rowMeans(tasic_trimmed[ , which(tasic_cell_types %in% list_other_cells[[i]])])
		if(i == 1){
			mean_others = data.frame(mean_other)
		} else {
			mean_others = cbind(mean_others, mean_other)
		}
    str(mean_others)
	}
  mean_fc = rowMeans(tasic_trimmed[ , which(tasic_cell_types %in% cell_type)]) /
      rowMeans(mean_others)
	mean_fc_names = data.frame(row.names(tasic_trimmed), mean_fc)
	names(mean_fc_names) = c("tasic_names", "mean_fc")

	#since already sorted, need to merge via rownames
	toptable = merge(tT, mean_fc_names, by.x = "row.names", by.y = "tasic_names")
	# toptable = toptable[toptable$mean_fc > fc_thres, ]
	toptable = toptable[order(toptable$t, decreasing = TRUE), ]

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

#create a modified volcano plot
plot(log(tasic_oli_df$mean_fc, 2), -log(tasic_oli_df$P.Value, 10))
plot(log(tasic_neu_df$mean_fc, 2), -log(tasic_neu_df$P.Value, 10))
plot(log(tasic_ast_df$mean_fc, 2), -log(tasic_ast_df$P.Value, 10))
plot(log(tasic_mic_df$mean_fc, 2), -log(tasic_mic_df$P.Value, 10))
plot(log(tasic_end_df$mean_fc, 2), -log(tasic_end_df$P.Value, 10))
