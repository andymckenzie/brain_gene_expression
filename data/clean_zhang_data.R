library(readxl)
library(zoo)
library(limma)
library(DGCA)

cutoff_thresh_percentile = 0.75

setwd("/Users/amckenz/Documents/github/brain_gene_expression/data/")

cell = read_excel("TableS4-HumanMouseMasterFPKMList.xlsx", col_names = FALSE, sheet = 2)

#http://stackoverflow.com/a/7735681/560791
#carry the last non-NA value forward
celltypes = unlist(cell[1, ])
celltypes = na.locf(celltypes[-1])
cellids = cell[2, ]
gender = cell[3, ]

combo = paste(celltypes, gender[-1], cellids[-1], sep = "_")
#remove identifier rows, gene column, and NA trailing columns
gnxp = cell[-c(1:3, 23227, 23228), -1]

gnxp = apply(gnxp, 2, as.numeric)
zhang_gnxp = gnxp
rownames(zhang_gnxp) = make.unique(cell[-c(1:3, 23227, 23228), 1])
colnames(zhang_gnxp) = celltypes
# pretty reasonably normalized
# plot(colMeans(as.matrix(zhang_gnxp_log)), colMedians(as.matrix(zhang_gnxp_log)))

zhang_gnxp_astr = zhang_gnxp[ , colnames(zhang_gnxp) %in% c("Human mature astrocytes")]
zhang_gnxp_neur = zhang_gnxp[ , colnames(zhang_gnxp) %in% c("Human Neurons")]
zhang_gnxp_olig = zhang_gnxp[ , colnames(zhang_gnxp) %in% c("Human Oligodendrocytes")]
zhang_gnxp_micr = zhang_gnxp[ , colnames(zhang_gnxp) %in% c("Human Microglia/Macrophage")]
zhang_gnxp_endo = zhang_gnxp[ , colnames(zhang_gnxp) %in% c("Human Endothelial")]

celltypes = gsub("Human mature astrocytes", "AST", celltypes)
celltypes = gsub("Human Neurons", "NEU", celltypes)
celltypes = gsub("Human Oligodendrocytes", "OLI", celltypes)
celltypes = gsub("Human Microglia/Macrophage", "MIC", celltypes)
celltypes = gsub("Human Endothelial", "END", celltypes)
celltypes = make.names(celltypes)

zhang_gnxp_log = log(zhang_gnxp + 1, 2)
zhang_gnxp_log = as.data.frame(zhang_gnxp_log)
fl = as.factor(celltypes)

dfxp_function <- function(cell_type, list_other_cells){

	#for a given cell type that you are considering
	#remove that gene if it has too low of expression *in that cell type*
	#find the average of each gene in its own cell type
  if(!cell_type == "NEU"){
    expr_average = rowMeans(zhang_gnxp_log[ , which(celltypes %in%
  		cell_type)])
  } else {
    expr_average = zhang_gnxp_log[ , which(celltypes %in% cell_type)]
  }
	cutoff = as.numeric(quantile(expr_average, cutoff_thresh_percentile, names = FALSE))
	zhang_trimmed = zhang_gnxp_log[(expr_average > cutoff), ]

  celltypes_contrast = celltypes
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
  print(design)
	fit = lmFit(zhang_trimmed, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
  str(fit2)
	fit2 = eBayes(fit2, 0.01, trend = TRUE)
  print(head(fit2$coefficients))
	tT = topTable(fit2, adjust = "BH", sort.by="B", number = nrow(fit2), coef = 1)

  #calculate the average vs all of the other cell types
	for(i in 1:length(list_other_cells)){
    #because the neuron cell type only contains one sample
    if(list_other_cells[[i]] == "NEU"){
      mean_other = zhang_trimmed[ , which(celltypes %in% list_other_cells[[i]])]
    } else {
      mean_other = rowMeans(zhang_trimmed[ , which(celltypes %in% list_other_cells[[i]])])
    }
		if(i == 1){
			mean_others = data.frame(mean_other)
		} else {
			mean_others = cbind(mean_others, mean_other)
		}
    str(mean_others)
	}
  if(!cell_type == "NEU"){
    mean_fc = rowMeans(zhang_trimmed[ , which(celltypes %in% cell_type)]) /
      rowMeans(mean_others)
  } else {
    mean_fc = zhang_trimmed[ , which(celltypes %in% cell_type)] /
      rowMeans(mean_others)
  }
	mean_fc_names = data.frame(row.names(zhang_trimmed), mean_fc)
	names(mean_fc_names) = c("zhang_names", "mean_fc")

	#since already sorted, need to merge via rownames
	toptable = merge(tT, mean_fc_names, by.x = "row.names", by.y = "zhang_names")
	toptable = toptable[order(toptable$P.Value), ]

	return(toptable)

}

######################################
# compare the main cell types to one another

z16_oli_df = dfxp_function("OLI", list("NEU", "AST", "MIC", "END"))
z16_neu_df = dfxp_function("NEU", list("OLI", "AST", "MIC", "END"))
z16_ast_df = dfxp_function("AST", list("NEU", "OLI", "MIC", "END"))
z16_mic_df = dfxp_function("MIC", list("NEU", "AST", "OLI", "END"))
z16_end_df = dfxp_function("END", list("NEU", "AST", "MIC", "OLI"))

#create a modified volcano plot
plot(z16_oli_df$mean_fc, -log(z16_oli_df$P.Value, 10))
plot(z16_neu_df$mean_fc, -log(z16_neu_df$P.Value, 10))
plot(z16_ast_df$mean_fc, -log(z16_ast_df$P.Value, 10))
plot(z16_mic_df$mean_fc, -log(z16_mic_df$P.Value, 10))
plot(z16_end_df$mean_fc, -log(z16_end_df$P.Value, 10))
