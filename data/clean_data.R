library(readxl)
library(edgeR)

############################
# functions

std_error <- function(x) sd(x)/sqrt(length(x))

#https://www.biostars.org/p/72846/
cpm_tmm <- function(counts){
    d = DGEList(counts = counts)
    d = calcNormFactors(d, method = "TMM")
    return(cpm(d, normalized.lib.sizes = TRUE))
}

################################
# sharma
#read in the sharma RPKM data set
sharma = read_excel("sharma_rna_nn.4160-S5.xlsx")
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

saveRDS(sharma_num, "data/sharma_rpkm.rds")

#read in the sharma proteomics data
sharma_prot = read_excel("data/sharma_proteomics_nn.4160-S7.xlsx", skip = 1)
sharma_prot_num = sharma_prot[, !colnames(sharma_prot) %in%
  c("Gene names", "Protein names", "PEP", "Mol. weight [kDa]",
  "Sequence coverage [%]", "Protein IDs", "Majority protein IDs")]
sharma_prot$`Gene names`[is.na(sharma_prot$`Gene names`)] = ""
rownames(sharma_prot_num) = make.unique(sharma_prot$`Gene names`)
#proteomics data is also approximately normalized... colMeans(as.matrix(sharma_num))

saveRDS(sharma_prot_num, "data/sharma_lfq.rds")

###################################
# darmanis

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

if(!identical(colnames(dar_norm_cells), type_gsm_clean_names)) stop("Pheno table doesn't match sample names.")

find_cell_mean_darmanis <- function(cell_type){

  tmp_gnxp = dar_norm_cells[ , dar_cell_types == cell_type]
  tmp_gnxp = as.matrix(tmp_gnxp)
  mean_gnxp = rowMeans(tmp_gnxp)
  se_gnxp = apply(tmp_gnxp, 1, std_error)

  tmp_gnxp_log = log(tmp_gnxp + 1, 2)
  mean_gnxp_log = rowMeans(tmp_gnxp_log)
  se_gnxp_log = apply(tmp_gnxp_log, 1, std_error)

  tmp_df = data.frame(mean_gnxp, se_gnxp, mean_gnxp_log, se_gnxp_log)

  colnames(tmp_df) = c(paste0(cell_type, "_mean"), paste0(cell_type, "_se"),
    paste0(cell_type, "_mean_log"), paste0(cell_type, "_se_log"))

  return(tmp_df)
}

oligodendrocytes_df = find_cell_mean_darmanis("Oligodendrocyte")
astrocytes_df = find_cell_mean_darmanis("Astrocyte")
endothelial_df = find_cell_mean_darmanis("Endothelial")
microglia_df = find_cell_mean_darmanis("Microglia")
neuron_df = find_cell_mean_darmanis("Neuron")
OPC_df = find_cell_mean_darmanis("OPC")
quiescent_fetal_df = find_cell_mean_darmanis("Quiescent Fetal Neuron")
replicating_fetal_df = find_cell_mean_darmanis("Replicating Fetal Neuron")

total_df = cbind(oligodendrocytes_df, astrocytes_df, endothelial_df,
  microglia_df, neuron_df, OPC_df, quiescent_fetal_df, replicating_fetal_df)
rownames(total_df) = rownames(dar_norm_cells)

saveRDS(total_df, "data/darmanis.rds")

###################################
# marques

gnxp = read.delim("data/GSE75330_Marques_et_al_mol_counts2.tab")
rownames(gnxp) = gnxp[ ,1]
gnxp = gnxp[ ,-1]
marques_norm = cpm_tmm(gnxp)

metadata = read.delim("data/metadate_oligos_in_tsne_22-Feb-2016.txt")

gnxp_cols = colnames(marques_norm)
cols_match_start = sapply(strsplit(gnxp_cols, ".", fixed = TRUE), "[", c(2,3))
cols_match_end = sapply(strsplit(gnxp_cols, ".", fixed = TRUE), "[", c(4))
cols_match_res = paste0(cols_match_start[1,], cols_match_start[2,], "_", cols_match_end)

metadata_match = metadata[match(cols_match_res, metadata$cellID), ]

stopifnot(identical(cols_match_res, metadata_match$cellID))

metadata_match$lev2_group = gsub("PPR", "VLMC", metadata_match$lev2_group)

find_cell_mean_marques <- function(cell_type){

  tmp_gnxp = marques_norm[ , metadata_match$lev2_group == cell_type]
  str(tmp_gnxp)
  tmp_gnxp = as.matrix(tmp_gnxp)
  mean_gnxp = rowMeans(tmp_gnxp)
  str(mean_gnxp)
  se_gnxp = apply(tmp_gnxp, 1, std_error)

  tmp_gnxp_log = log(tmp_gnxp + 1, 2)
  mean_gnxp_log = rowMeans(tmp_gnxp_log)
  se_gnxp_log = apply(tmp_gnxp_log, 1, std_error)

  tmp_df = data.frame(mean_gnxp, se_gnxp, mean_gnxp_log, se_gnxp_log)

  colnames(tmp_df) = c(paste0(cell_type, "_mean"), paste0(cell_type, "_se"),
    paste0(cell_type, "_mean_log"), paste0(cell_type, "_se_log"))

  return(tmp_df)
}

oligo_types = unique(metadata_match$lev2_group)
for(i in 1:length(oligo_types)){
  tmp_df = find_cell_mean_marques(oligo_types[i])
  if(i == 1){
    marques_total_df = tmp_df
  } else {
    marques_total_df = cbind(marques_total_df, tmp_df)
  }
  print(i)
}

rownames(marques_total_df) = rownames(marques_norm)
saveRDS(marques_total_df, "data/marques.rds")

###############################
# tasic

metadata = read.csv("data/tasic_cell_metadata.csv")
counts = read.csv("data/tasic_genes_counts.csv")
rownames(counts) = counts[ , 1]
counts = counts[ , -1]
tasic_norm = cpm_tmm(counts)

find_cell_mean_tasic <- function(cell_type){

  tmp_gnxp = tasic_norm[ , metadata$sub_class == cell_type]
  tmp_gnxp = as.matrix(tmp_gnxp)
  mean_gnxp = rowMeans(tmp_gnxp)
  str(mean_gnxp)
  se_gnxp = apply(tmp_gnxp, 1, std_error)

  tmp_gnxp_log = log(tmp_gnxp + 1, 2)
  mean_gnxp_log = rowMeans(tmp_gnxp_log)
  se_gnxp_log = apply(tmp_gnxp_log, 1, std_error)

  tmp_df = data.frame(mean_gnxp, se_gnxp, mean_gnxp_log, se_gnxp_log)

  colnames(tmp_df) = c(paste0(cell_type, "_mean"), paste0(cell_type, "_se"),
    paste0(cell_type, "_mean_log"), paste0(cell_type, "_se_log"))

  return(tmp_df)
}

sub_class = table(metadata$sub_class)
sub_class = sub_class[sub_class >= 15]
sub_class = sub_class[!names(sub_class) == ""]
cell_types = names(sub_class)

for(i in 1:length(cell_types)){
  tmp_df = find_cell_mean_tasic(cell_types[i])
  if(i == 1){
    tasic_total_df = tmp_df
  } else {
    tasic_total_df = cbind(tasic_total_df, tmp_df)
  }
  print(i)
}

colnames(tasic_total_df)[grepl("Sst", colnames(tasic_total_df))] =
  c("Sst Cbln4_mean", "Sst Cbln4_se", "Sst Cbln4_mean_log", "Sst Cbln4_se_log",
  "Sst Chodl_mean", "Sst Chodl_se", "Sst Chodl_mean_log", "Sst Chodl_se_log",
  "Sst Nr2f2_mean", "Sst Nr2f2_se", "Sst Nr2f2_mean_log", "Sst Nr2f2_se_log")

rownames(tasic_total_df) = rownames(tasic_norm)
saveRDS(tasic_total_df, "data/tasic.rds")

################################
# zeisel

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
cell_types = gsub("Mgl2", "Microglia", cell_types)
cell_types = gsub("Mgl2", "Microglia", cell_types)
cell_types[grepl("Int", cell_types)] = "Interneuron"
cell_types[grepl("Oligo", cell_types)] = "Oligodendrocyte"
cell_types = gsub("Vsmc", "Vascular Smooth Muscle", cell_types)
cell_types = gsub("Peric", "Pericyte", cell_types)
cell_types = gsub("Epend", "Ependymal", cell_types)
cell_types = gsub("Pyr", " Pyramidal ", cell_types)
cell_types = gsub("Sub Pyramidal ", "Subiculum Pyramidal", cell_types)
cell_types = gsub("Pyramidal DL", "Pyramidal Deep Layer", cell_types)

cell_types_table = table(cell_types)
cell_types_include = cell_types_table[cell_types_table >= 15]
cell_types_include = cell_types_include[!names(cell_types_include) == "(none)"]

find_cell_mean_zeisel <- function(cell_type){

  tmp_gnxp = norm_zeis[ , cell_types == cell_type]
  tmp_gnxp = as.matrix(tmp_gnxp)
  mean_gnxp = rowMeans(tmp_gnxp)
  str(mean_gnxp)
  se_gnxp = apply(tmp_gnxp, 1, std_error)

  tmp_gnxp_log = log(tmp_gnxp + 1, 2)
  mean_gnxp_log = rowMeans(tmp_gnxp_log)
  se_gnxp_log = apply(tmp_gnxp_log, 1, std_error)

  tmp_df = data.frame(mean_gnxp, se_gnxp, mean_gnxp_log, se_gnxp_log)

  colnames(tmp_df) = c(paste0(cell_type, "_mean"), paste0(cell_type, "_se"),
    paste0(cell_type, "_mean_log"), paste0(cell_type, "_se_log"))

  return(tmp_df)
}

for(i in 1:length(cell_types_include)){
  tmp_df = find_cell_mean_zeisel(names(cell_types_include)[i])
  if(i == 1){
    zeisel_total_df = tmp_df
  } else {
    zeisel_total_df = cbind(zeisel_total_df, tmp_df)
  }
  print(i)
}

rownames(zeisel_total_df) = rownames(norm_zeis)
saveRDS(zeisel_total_df, "data/zeisel.rds")
