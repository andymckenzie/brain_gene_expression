library(readxl)
library(zoo)
library(limma)
library(DGCA)
library(HGNChelper)


setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

cell = read_excel("data/TableS4-HumanMouseMasterFPKMList.xlsx", col_names = FALSE, sheet = 2)

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

##########################
#convert to cleaned human gene sybmols where possible
source("/Users/amckenz/Documents/github/alzolig/convert_mgi_to_hgnc.R")
gene_symbols = switchGenesToHGCN(rownames(zhang_gnxp_log))
rownames(zhang_gnxp_log) = make.unique(gene_symbols)
collapsed_df = collapseRows(zhang_gnxp_log, rowGroup = gene_symbols, rowID = make.unique(gene_symbols), method="MaxMean")
zhang_gnxp_log = zhang_gnxp_log[collapsed_df$selectedRow, ]
gene_symbols = gene_symbols[collapsed_df$selectedRow]
rownames(zhang_gnxp_log) = gene_symbols

dfxp_function <- function(cell_type, list_other_cells){

  if(!cell_type == "NEU"){
    expr_average = rowMeans(zhang_gnxp_log[ , which(celltypes %in% cell_type)])
  } else {
    expr_average = zhang_gnxp_log[ , which(celltypes %in% cell_type)]
  }

  celltypes_contrast = celltypes
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
	fit = lmFit(zhang_gnxp_log, design)
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

######################################
# compare the main cell types to one another

z16_oli_df = dfxp_function("OLI", list("NEU", "AST", "MIC", "END"))
z16_neu_df = dfxp_function("NEU", list("OLI", "AST", "MIC", "END"))
z16_ast_df = dfxp_function("AST", list("NEU", "OLI", "MIC", "END"))
z16_mic_df = dfxp_function("MIC", list("NEU", "AST", "OLI", "END"))
z16_end_df = dfxp_function("END", list("NEU", "AST", "MIC", "OLI"))

saveRDS(z16_oli_df, "data/dfxp/z16_oli_df_dfxp.rds")
saveRDS(z16_neu_df, "data/dfxp/z16_neu_df_dfxp.rds")
saveRDS(z16_ast_df, "data/dfxp/z16_ast_df_dfxp.rds")
saveRDS(z16_mic_df, "data/dfxp/z16_mic_df_dfxp.rds")
saveRDS(z16_end_df, "data/dfxp/z16_end_df_dfxp.rds")
