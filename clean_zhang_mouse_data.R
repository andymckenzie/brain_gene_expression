#oops this is the same as the other data set 

library(readxl)
library(zoo)
library(limma)
library(DGCA)

cutoff_thresh_percentile = 0.75

setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

mouse_cell = read_excel("TableS4-HumanMouseMasterFPKMList.xlsx", col_names = FALSE, sheet = 3)

#http://stackoverflow.com/a/7735681/560791
#carry the last non-NA value forward
celltypes = unlist(mouse_cell[1, ])
celltypes = na.locf(celltypes[-1])
#remove identifier rows, gene column
gnxp = mouse_cell[-c(1:2), -1]

gnxp = apply(gnxp, 2, as.numeric)
zhang_gnxp = gnxp
rownames(zhang_gnxp) = make.unique(mouse_cell[-c(1:2), 1])
colnames(zhang_gnxp) = celltypes
zhang_gnxp_log = log(zhang_gnxp + 1, 2)
#seems reasonably normalized
# boxplot(as.matrix(zhang_gnxp_log))

#pooling across ages for AST
#1 month	4 months	7 months	9 months
celltypes = gsub("Mouse Hepacam immunopanned astrocytes", "AST", celltypes)
#Neuron 3	Neuron 4
celltypes = gsub("Mouse neurons", "NEU", celltypes)
#Oligodendrocyte precursor cell 3	Oligodendrocyte precursor cell 4	Newly formed oligodendrocyte 3	Newly formed oligodendrocyte 4	Myelinating oligodendrocyte 4	Myelinating oligodenrocyte 5
oligos = grep("oligodendrocytes", celltypes)
celltypes[oligos[1:2]] = "OPC"
celltypes[oligos[3:4]] = "NFO"
celltypes[oligos[5:6]] = "MOL"
celltypes = gsub("Mouse microglia", "MIC", celltypes)
celltypes = gsub("Mouse endothelial", "END", celltypes)
celltypes = make.names(celltypes)

dfxp_function <- function(cell_type, list_other_cells){

	#for a given cell type that you are considering
	#remove that gene if it has too low of expression *in that cell type*
	#find the average of each gene in its own cell type
  expr_average = rowMeans(zhang_gnxp_log[ , which(celltypes %in% cell_type)])
  print(sum(is.na(expr_average)))
	cutoff = as.numeric(quantile(expr_average, cutoff_thresh_percentile, names = FALSE))
  str(cutoff)
	zhang_trimmed = zhang_gnxp_log[(expr_average > cutoff), ]
  str(zhang_trimmed)

  celltypes_contrast = celltypes
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
  str(design)
  str(zhang_trimmed)
	fit = lmFit(zhang_trimmed, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
  str(fit2)
	fit2 = eBayes(fit2, 0.01, trend = TRUE)
  print(head(fit2$coefficients))
	tT = topTable(fit2, adjust = "BH", sort.by="B", number = nrow(fit2), coef = 1, confint = TRUE)

  #calculate the average vs all of the other cell types
	for(i in 1:length(list_other_cells)){
    mean_other = rowMeans(zhang_trimmed[ , which(celltypes %in% list_other_cells[[i]])])
		if(i == 1){
			mean_others = data.frame(mean_other)
		} else {
			mean_others = cbind(mean_others, mean_other)
		}
    str(mean_others)
	}
  mean_fc = rowMeans(zhang_trimmed[ , which(celltypes %in% cell_type)]) / rowMeans(mean_others)
	mean_fc_names = data.frame(row.names(zhang_trimmed), mean_fc)
	names(mean_fc_names) = c("zhang_names", "mean_fc")

	#since already sorted, need to merge via rownames
  str(tT)
  str(mean_fc_names)
	toptable = merge(tT, mean_fc_names, by.x = "row.names", by.y = "zhang_names")
	toptable = toptable[order(toptable$t, decreasing = TRUE), ]

	return(toptable)

}

######################################
# compare the main cell types to one another

zm16_oli_df = dfxp_function("MOL", list("NEU", "AST", "MIC", "END"))
zm16_neu_df = dfxp_function("NEU", list("MOL", "AST", "MIC", "END", "OPC"))
zm16_ast_df = dfxp_function("AST", list("NEU", "MOL", "MIC", "END", "OPC"))
zm16_mic_df = dfxp_function("MIC", list("NEU", "AST", "MOL", "END", "OPC"))
zm16_end_df = dfxp_function("END", list("NEU", "AST", "MIC", "MOL", "OPC"))
zm16_opc_df = dfxp_function("OPC", list("NEU", "AST", "MIC", "END"))

saveRDS(zm16_oli_df, "data/dfxp/zm16_oli_df_dfxp.rds")
saveRDS(zm16_neu_df, "data/dfxp/zm16_neu_df_dfxp.rds")
saveRDS(zm16_ast_df, "data/dfxp/zm16_ast_df_dfxp.rds")
saveRDS(zm16_mic_df, "data/dfxp/zm16_mic_df_dfxp.rds")
saveRDS(zm16_end_df, "data/dfxp/zm16_end_df_dfxp.rds")
saveRDS(zm16_opc_df, "data/dfxp/zm16_opc_df_dfxp.rds")
