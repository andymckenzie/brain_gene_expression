
library(edgeR)
library(DGCA)
library(limma)
library(HGNChelper)


############################
# functions

#https://www.biostars.org/p/72846/
cpm_tmm <- function(counts){
    d = DGEList(counts = counts)
    d = calcNormFactors(d, method = "TMM")
    return(cpm(d, normalized.lib.sizes = TRUE))
}

###################################
# darmanis

setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

filenames = list.files("data/GSE67835_RAW",
  pattern = "*.csv.gz", full.names = TRUE)

#may get Error: OutOfMemoryError (Java): Java heap space due to XLConnect here
ldf = lapply(filenames, read.delim, header = FALSE)
#this function takes time -- on the scale of minutes
dar_merged = Reduce(function(...) merge(..., by = "V1", all = TRUE),
	ldf)

sample_names = list.files("data/GSE67835_RAW", pattern="*.csv.gz")
sample_names = sapply(strsplit(sample_names, "_"), "[[", 1)

rownames(dar_merged) = dar_merged[ , 1]
dar_merged = dar_merged[ , -1]
colnames(dar_merged) = sample_names
rownames(dar_merged) = gsub(" ", "", rownames(dar_merged))
dar_norm = cpm_tmm(dar_merged)

samples = read.delim("data/GSE67835-GPL18573_series_matrix.txt", fill = TRUE, skip = 38, quote = "", header = T)
samples2 = read.delim("data/GSE67835-GPL15520_series_matrix.txt", fill = TRUE, skip = 38, quote = "", header = T)
samples_full = cbind(samples, samples2)

#find the GSM's that correspond to each of the cell types
types = samples_full[9, ]
type_gsm = as.character(names(types))
type_gsm_clean = sapply(strsplit(type_gsm, ".", fixed = TRUE), "[[", 2)
type_gsm_clean_names = type_gsm_clean[!grepl("Sample_geo_accession", type_gsm_clean)]
dar_cell_types = as.character(types)[!grepl("Sample_geo_accession", type_gsm_clean)]
dar_norm_cells = dar_norm[ , match(type_gsm_clean_names, colnames(dar_norm))]

dar_cell_types = gsub("\"cell type: oligodendrocytes\"", "Oligodendrocyte", dar_cell_types)
dar_cell_types = gsub("\"cell type: astrocytes\"", "Astrocyte", dar_cell_types)
dar_cell_types = gsub("\"cell type: endothelial\"", "Endothelial", dar_cell_types)
dar_cell_types = gsub("\"cell type: microglia\"", "Microglia", dar_cell_types)
dar_cell_types = gsub("\"cell type: neurons\"", "Neuron", dar_cell_types)
dar_cell_types = gsub("\"cell type: OPC\"", "OPC", dar_cell_types)
dar_cell_types = gsub("\"cell type: fetal_quiescent\"", "Quiescent Fetal Neuron", dar_cell_types)
dar_cell_types = gsub("\"cell type: fetal_replicating\"", "Replicating Fetal Neuron", dar_cell_types)
dar_cell_types = make.names(dar_cell_types)

if(!identical(colnames(dar_norm_cells), type_gsm_clean_names)) stop("Pheno table doesn't match sample names.")

dar_gnxp_log = log(dar_norm_cells + 1, 2)

##########################
#convert to cleaned human gene sybmols where possible
source("/Users/amckenz/Documents/github/alzolig/convert_mgi_to_hgnc.R")
gene_symbols = switchGenesToHGCN(rownames(dar_merged))
rownames(dar_merged) = make.unique(gene_symbols)
collapsed_df = collapseRows(dar_merged, rowGroup = gene_symbols, rowID = make.unique(gene_symbols), method="MaxMean")
dar_merged = dar_merged[collapsed_df$selectedRow, ]
dar_norm_cells = dar_norm_cells[collapsed_df$selectedRow, ]
gene_symbols = gene_symbols[collapsed_df$selectedRow]
rownames(dar_merged) = gene_symbols

#############################
# dfxp

dfxp_function <- function(cell_type, list_other_cells){

  expr_average = rowMeans(as.matrix(dar_norm_cells[ , dar_cell_types == cell_type]))

  celltypes_contrast = dar_cell_types
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)

  dar_voom_res = voom(dar_merged, design, plot = TRUE)
  fit = lmFit(dar_voom_res, design)
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

d15_oli_df = dfxp_function("Oligodendrocyte", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))
d15_neu_df = dfxp_function("Neuron", list("Oligodendrocyte", "Astrocyte", "Microglia", "Endothelial", "OPC"))
d15_ast_df = dfxp_function("Astrocyte", list("Neuron", "Oligodendrocyte", "Microglia", "Endothelial", "OPC"))
d15_mic_df = dfxp_function("Microglia", list("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial", "OPC"))
d15_end_df = dfxp_function("Endothelial", list("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte", "OPC"))
d15_opc_df = dfxp_function("OPC", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))

saveRDS(d15_oli_df, "data/dfxp/d15_oli_df_dfxp.rds")
saveRDS(d15_neu_df, "data/dfxp/d15_neu_df_dfxp.rds")
saveRDS(d15_ast_df, "data/dfxp/d15_ast_df_dfxp.rds")
saveRDS(d15_mic_df, "data/dfxp/d15_mic_df_dfxp.rds")
saveRDS(d15_end_df, "data/dfxp/d15_end_df_dfxp.rds")
saveRDS(d15_opc_df, "data/dfxp/d15_opc_df_dfxp.rds")
