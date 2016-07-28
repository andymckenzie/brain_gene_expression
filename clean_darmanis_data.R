
library(edgeR)
library(DGCA)

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
# darmanis

cutoff_thresh_percentile = 0.5

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

dfxp_function <- function(cell_type, list_other_cells){

	#for a given cell type that you are considering
	#remove that gene if it has too low of expression *in that cell type*
	#find the average of each gene in its own cell type
  expr_average = rowMeans(as.matrix(dar_norm_cells[ , dar_cell_types == cell_type]))
	cutoff = as.numeric(quantile(expr_average, cutoff_thresh_percentile, names = FALSE))
	dar_trimmed = dar_gnxp_log[(expr_average > cutoff), ]

  celltypes_contrast = dar_cell_types
  celltypes_contrast = gsub(cell_type, "MAIN", celltypes_contrast)
  for(i in 1:length(list_other_cells)){
    celltypes_contrast = gsub(list_other_cells[[i]], "OTHERS", celltypes_contrast)
  }
  design = makeDesign(celltypes_contrast)
	fit = lmFit(dar_trimmed, design)
  contrast_matrix = makeContrasts(MAIN-OTHERS, levels = as.factor(celltypes_contrast))

	fit2 = contrasts.fit(fit, contrast_matrix)
	fit2 = eBayes(fit2, 0.01, trend = TRUE)
  print(head(fit2$coefficients))
	tT = topTable(fit2, adjust = "BH", sort.by="B", number = nrow(fit2), coef = 1)
	# tT = tT[tT$P.Value < pval_thresh, ]

  #calculate the average vs all of the other cell types
	for(i in 1:length(list_other_cells)){
    mean_other = rowMeans(dar_trimmed[ , which(dar_cell_types %in% list_other_cells[[i]])])
		if(i == 1){
			mean_others = data.frame(mean_other)
		} else {
			mean_others = cbind(mean_others, mean_other)
		}
    str(mean_others)
	}
  mean_fc = rowMeans(dar_trimmed[ , which(dar_cell_types %in% cell_type)]) /
      rowMeans(mean_others)
	mean_fc_names = data.frame(row.names(dar_trimmed), mean_fc)
	names(mean_fc_names) = c("dar_names", "mean_fc")

	#since already sorted, need to merge via rownames
	toptable = merge(tT, mean_fc_names, by.x = "row.names", by.y = "dar_names")
	toptable = toptable[order(toptable$t, decreasing = TRUE), ]

	return(toptable)

}

d15_oli_df = dfxp_function("Oligodendrocyte", list("Neuron", "Astrocyte", "Microglia", "Endothelial"))
d15_neu_df = dfxp_function("Neuron", list("Oligodendrocyte", "Astrocyte", "Microglia", "Endothelial"))
d15_ast_df = dfxp_function("Astrocyte", list("Neuron", "Oligodendrocyte", "Microglia", "Endothelial"))
d15_mic_df = dfxp_function("Microglia", list("Neuron", "Astrocyte", "Oligodendrocyte", "Endothelial"))
d15_end_df = dfxp_function("Endothelial", list("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"))

#create a modified volcano plot
plot(log(d15_oli_df$mean_fc, 2), -log(d15_oli_df$P.Value, 10))
plot(log(d15_neu_df$mean_fc, 2), -log(d15_neu_df$P.Value, 10))
plot(log(d15_ast_df$mean_fc, 2), -log(d15_ast_df$P.Value, 10))
plot(log(d15_mic_df$mean_fc, 2), -log(d15_mic_df$P.Value, 10))
plot(log(d15_end_df$mean_fc, 2), -log(d15_end_df$P.Value, 10))
