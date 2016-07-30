
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
# zeis

zeis = read.delim("data/ziesel_expression_mRNA_17-Aug-2014.txt")
zeis_gnxp = zeis[-c(1:10), ]
rownames(zeis_gnxp) = make.unique(zeis_gnxp[ , 1])
zeis_gnxp = zeis_gnxp[ , -c(1:2)]
zeis_gnxp = data.matrix(zeis_gnxp)
norm_zeis = cpm_tmm(zeis_gnxp)

cell_types = as.character(zeis[9,])[-c(1,2)]
cell_types = gsub("Astro1", "Astrocyte", cell_types)
cell_types = gsub("Astro2", "Astrocyte", cell_types)
cell_types = gsub("Vend1", "Endothelial", cell_types)
cell_types = gsub("Vend2", "Endothelial", cell_types)
cell_types = gsub("Pvm1", "Perivascular macrophage", cell_types)
cell_types = gsub("Pvm2", "Perivascular macrophage", cell_types)
cell_types = gsub("Mgl1", "Microglia", cell_types)
cell_types = gsub("Mgl2", "Microglia", cell_types)
cell_types[grepl("Int", cell_types)] = "Interneuron"
#likely representing stages of maturation:
# immature (Oligo4), pre-myelinating (Oligo2), myelinating
# (Oligo5) and terminally differentiated post-myelination
# (Oligo6) oligodendrocytes. An intermediate population, Oligo3,
# was almost exclusively observed in somatosensory cortex
# and may represent a distinct cellular state specific for
# this tissue. The subclass Oligo1, which did not express the
# prototypical genes associated with oligodendrocyte precursor
# cells (OPCs), may represent a postmitotic cellular state,
# associated with the first steps of oligodendrocyte differentiation.
# Together, the Oligo1 – Oligo6
# populations may represent sequential steps in the process
# of maturation from an OPC to a terminally
# differentiated oligodendrocyte.
cell_types = gsub("Oligo4", "OPC", cell_types)
cell_types = gsub("Oligo5", "Oligodendrocyte", cell_types)
cell_types = gsub("Oligo6", "Oligodendrocyte", cell_types)
cell_types = gsub("Vsmc", "Vascular Smooth Muscle", cell_types)
cell_types = gsub("Peric", "Pericyte", cell_types)
cell_types = gsub("Epend", "Ependymal", cell_types)
cell_types = gsub("Pyr", "Pyramidal", cell_types)
cell_types = gsub("Sub Pyramidal ", "Subiculum Pyramidal", cell_types)
cell_types = gsub("Pyramidal DL", "Pyramidal Deep Layer", cell_types)

cell_types_main = cell_types
cell_types_main = gsub("Interneuron", "Neuron", cell_types_main)
cell_types_main[grepl("Pyramidal", cell_types_main)] =  "Neuron"
cell_types_main = make.names(cell_types_main)

norm_zeis_log = log(norm_zeis + 1, 2)

dfxp_function <- function(cell_type, list_other_cells){

	#for a given cell type that you are considering
	#remove that gene if it has too low of expression *in that cell type*
	#find the average of each gene in its own cell type
  expr_average = rowMeans(as.matrix(norm_zeis_log[ , cell_types_main == cell_type]))
	cutoff = as.numeric(quantile(expr_average, cutoff_thresh_percentile, names = FALSE))
	zeis_trimmed = norm_zeis_log[(expr_average > cutoff), ]

  #The voom transformation
# is applied to the read counts. This converts the counts to log-counts per million with associated
# precision weights. After this, the RNA-seq data can be analyzed as if it was microarray data

  celltypes_contrast = cell_types_main
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
	fit = lmFit(zeis_trimmed, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
	fit2 = eBayes(fit2, 0.01, trend = TRUE)
  print(head(fit2$coefficients))
	tT = topTable(fit2, adjust = "BH", sort.by="B", number = nrow(fit2), coef = 1, confint = TRUE)
	# tT = tT[tT$P.Value < pval_thresh, ]

  #calculate the average vs all of the other cell types
	for(i in 1:length(list_other_cells)){
    mean_other = rowMeans(zeis_trimmed[ , which(cell_types_main %in% list_other_cells[[i]])])
		if(i == 1){
			mean_others = data.frame(mean_other)
		} else {
			mean_others = cbind(mean_others, mean_other)
		}
    str(mean_others)
	}
  mean_fc = rowMeans(zeis_trimmed[ , which(cell_types_main %in% cell_type)]) /
      rowMeans(mean_others)
	mean_fc_names = data.frame(row.names(zeis_trimmed), mean_fc)
	names(mean_fc_names) = c("zeis_names", "mean_fc")

	#since already sorted, need to merge via rownames
	toptable = merge(tT, mean_fc_names, by.x = "row.names", by.y = "zeis_names")
	# toptable = toptable[toptable$mean_fc > fc_thres, ]
	toptable = toptable[order(toptable$t, decreasing = TRUE), ]

	return(toptable)

}

zeis_oli_df = dfxp_function("Oligodendrocyte", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))
zeis_neu_df = dfxp_function("Neuron", list("Oligodendrocyte", "Astrocyte", "Microglia", "Endothelial", "OPC"))
zeis_ast_df = dfxp_function("Astrocyte", list("Neuron", "Oligodendrocyte", "Microglia", "Endothelial", "OPC"))
zeis_mic_df = dfxp_function("Microglia", list("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial", "OPC"))
zeis_end_df = dfxp_function("Endothelial", list("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC"))
zeis_opc_df = dfxp_function("OPC", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))

saveRDS(zeis_oli_df, "data/dfxp/zeis_oli_df_dfxp.rds")
saveRDS(zeis_neu_df, "data/dfxp/zeis_neu_df_dfxp.rds")
saveRDS(zeis_ast_df, "data/dfxp/zeis_ast_df_dfxp.rds")
saveRDS(zeis_mic_df, "data/dfxp/zeis_mic_df_dfxp.rds")
saveRDS(zeis_end_df, "data/dfxp/zeis_end_df_dfxp.rds")
saveRDS(zeis_opc_df, "data/dfxp/zeis_opc_df_dfxp.rds")

#create a modified volcano plot
plot(log(zeis_oli_df$mean_fc, 2), -log(zeis_oli_df$P.Value, 10))
plot(log(zeis_neu_df$mean_fc, 2), -log(zeis_neu_df$P.Value, 10))
plot(log(zeis_ast_df$mean_fc, 2), -log(zeis_ast_df$P.Value, 10))
plot(log(zeis_mic_df$mean_fc, 2), -log(zeis_mic_df$P.Value, 10))
plot(log(zeis_end_df$mean_fc, 2), -log(zeis_end_df$P.Value, 10))
