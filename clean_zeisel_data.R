
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
# zeis

zeis = read.delim("data/ziesel_expression_mRNA_17-Aug-2014.txt")
zeis_gnxp = zeis[-c(1:10), ]
zeis_genes = zeis_gnxp[ , 1]
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
# cell_types = gsub("Oligo4", "OPC", cell_types)
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

##########################
#convert to cleaned human gene sybmols where possible
source("/Users/amckenz/Documents/github/alzolig/convert_mgi_to_hgnc.R")
gene_symbols = convert_mgi_to_hgnc(zeis_genes, genome = "mm10", return_df = FALSE)
gene_symbols = switchGenesToHGCN(gene_symbols)
rownames(zeis_gnxp) = make.unique(gene_symbols)
collapsed_df = collapseRows(zeis_gnxp, rowGroup = gene_symbols, rowID = make.unique(gene_symbols), method="MaxMean")
zeis_gnxp = zeis_gnxp[collapsed_df$selectedRow, ]
norm_zeis_log = norm_zeis_log[collapsed_df$selectedRow, ]
gene_symbols = gene_symbols[collapsed_df$selectedRow]
rownames(zeis_gnxp) = gene_symbols

####################
#dfxp

dfxp_function <- function(cell_type, list_other_cells){

  expr_average = rowMeans(as.matrix(norm_zeis_log[ , cell_types_main == cell_type]))

  celltypes_contrast = cell_types_main
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)

  #The voom transformation is applied to the read counts. This converts the counts to log-counts per million with associated
  # precision weights. After this, the RNA-seq data can be analyzed as if it was microarray data
  zeis_voom_res = voom(zeis_gnxp, design, plot = TRUE)
	fit = lmFit(zeis_voom_res, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
  #If, as you say, voom already handles this issue, excellent, there is no need to set trend=TRUE in eBayes function.
	fit2 = eBayes(fit2)
	toptable = topTable(fit2, adjust = "BH", sort.by = "none", number = nrow(fit2), coef = 1, confint = TRUE)

  toptable$expr_average = expr_average
  source("shrink_logFC.R")
  toptable$fc_zscore = shrink_pval_and_avg_expr(toptable)
  toptable$genes = rownames(toptable)
	toptable = toptable[order(toptable$fc_zscore, decreasing = TRUE), ]

	return(toptable)

}

zeis_oli_df = dfxp_function("Oligodendrocyte", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))
zeis_neu_df = dfxp_function("Neuron", list("Oligodendrocyte", "Astrocyte", "Microglia", "Endothelial"))
zeis_ast_df = dfxp_function("Astrocyte", list("Neuron", "Oligodendrocyte", "Microglia", "Endothelial"))
zeis_mic_df = dfxp_function("Microglia", list("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial"))
zeis_end_df = dfxp_function("Endothelial", list("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"))

saveRDS(zeis_oli_df, "data/dfxp/zeis_oli_df_dfxp.rds")
saveRDS(zeis_neu_df, "data/dfxp/zeis_neu_df_dfxp.rds")
saveRDS(zeis_ast_df, "data/dfxp/zeis_ast_df_dfxp.rds")
saveRDS(zeis_mic_df, "data/dfxp/zeis_mic_df_dfxp.rds")
saveRDS(zeis_end_df, "data/dfxp/zeis_end_df_dfxp.rds")
