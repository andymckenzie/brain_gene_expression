library(readxl)
library(edgeR)
library(DGCA)

cutoff_thresh_percentile = 0.5

setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

############################
# functions

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
sharma_cell_types[grepl("Oligodendrocytes div4", sharma_cell_types)] = "Oligodendrocyte"
sharma_cell_types[grepl("adult microglia", sharma_cell_types)] = "Microglia"
sharma_cell_types[grepl("cortical neurons div10", sharma_cell_types)] = "Neuron"
sharma_cell_types[grepl("Astrocytes", sharma_cell_types)] = "Astrocyte"
sharma_cell_types = make.names(sharma_cell_types)

sharma_num_log = log(sharma_num + 1, 2)

dfxp_function <- function(cell_type, list_other_cells){

	#for a given cell type that you are considering
	#remove that gene if it has too low of expression *in that cell type*
	#find the average of each gene in its own cell type
  expr_average = rowMeans(as.matrix(sharma_num_log[ , sharma_cell_types == cell_type]))
	cutoff = as.numeric(quantile(expr_average, cutoff_thresh_percentile, names = FALSE))
	sharma_trimmed = sharma_num_log[(expr_average > cutoff), ]

  celltypes_contrast = sharma_cell_types
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
	fit = lmFit(sharma_trimmed, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
	fit2 = eBayes(fit2, 0.01, trend = TRUE)
  print(head(fit2$coefficients))
	tT = topTable(fit2, adjust = "BH", sort.by="B", number = nrow(fit2), coef = 1)
	# tT = tT[tT$P.Value < pval_thresh, ]

  #calculate the average vs all of the other cell types
	for(i in 1:length(list_other_cells)){
    mean_other = rowMeans(sharma_trimmed[ , which(sharma_cell_types %in% list_other_cells[[i]])])
		if(i == 1){
			mean_others = data.frame(mean_other)
		} else {
			mean_others = cbind(mean_others, mean_other)
		}
    str(mean_others)
	}
  mean_fc = rowMeans(sharma_trimmed[ , which(sharma_cell_types %in% cell_type)]) /
      rowMeans(mean_others)
	mean_fc_names = data.frame(row.names(sharma_trimmed), mean_fc)
	names(mean_fc_names) = c("sharma_names", "mean_fc")

	#since already sorted, need to merge via rownames
	toptable = merge(tT, mean_fc_names, by.x = "row.names", by.y = "sharma_names")
	# toptable = toptable[toptable$mean_fc > fc_thres, ]
	toptable = toptable[order(toptable$t, decreasing = TRUE), ]

	return(toptable)

}

sharma_oli_df = dfxp_function("Oligodendrocyte", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))
sharma_neu_df = dfxp_function("Neuron", list("Oligodendrocyte", "Astrocyte", "Microglia", "Endothelial"))
sharma_ast_df = dfxp_function("Astrocyte", list("Neuron", "Oligodendrocyte", "Microglia", "Endothelial"))
sharma_mic_df = dfxp_function("Microglia", list("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial"))

#create a modified volcano plot
plot(log(sharma_oli_df$mean_fc, 2), -log(sharma_oli_df$P.Value, 10))
plot(log(sharma_neu_df$mean_fc, 2), -log(sharma_neu_df$P.Value, 10))
plot(log(sharma_ast_df$mean_fc, 2), -log(sharma_ast_df$P.Value, 10))
plot(log(sharma_mic_df$mean_fc, 2), -log(sharma_mic_df$P.Value, 10))
