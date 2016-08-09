
library(HGNChelper)
library(ggplot2)
library(DGCA)
library(bayesbio)
library(corrplot)
library(psych) #geometric.mean
library(gridExtra) #grid.arrange
library(gtable)

cbPalette = c("#CC79A7", "#E69F00", "#56B4E9", "#000000", "#009E73", "#F0E442", "#0072B2",
 "red", "lightgreen", "#999999", "#990000")

setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

source("/Users/amckenz/Dropbox/zhang/general_code/get_fold_enrichment.R")

#a function that takes two of the df's, orders and identifies top N with one of the cutoff types,
#then find jaccard index of the plot, and spit back the results.
#part of this function needs to clean the gene symbols with HGNChelper
#and then purge/merge the df's so that only the shared genes are included.
compare_two_signatures <- function(df1, df2, nTerms, orderByCol, orderDecr, cell){

  universe = intersect(df1$genes, df2$genes)
  df1 = df1[df1$genes %in% universe, ]
  df2 = df2[df2$genes %in% universe, ]
  df2 = df2[!duplicated(df2$genes), ]

  df1 = df1[order(df1[ , orderByCol], decreasing = orderDecr), ]
  df2 = df2[order(df2[ , orderByCol], decreasing = orderDecr), ]

  fe = vector()
  pvals = vector()

  for(i in 1:length(nTerms)){
    top_n_df1 = unique(head(df1, nTerms[i])$genes)
    top_n_df2 = unique(head(df2, nTerms[i])$genes)
    res = get_fold_enrichment(top_n_df1, top_n_df2, universe)
    fe[i] = res[[1]]
    pvals[i] = res[[2]]
  }

  return(data.frame(nTerms = nTerms, fe = fe, pvals = pvals, cell = rep(cell, length(nTerms))))

}

z15_oli_df = readRDS("data/dfxp/z15_oli_df_dfxp.rds")
z15_ast_df = readRDS("data/dfxp/z15_ast_df_dfxp.rds")
z15_neu_df = readRDS("data/dfxp/z15_neu_df_dfxp.rds")
z15_mic_df = readRDS("data/dfxp/z15_mic_df_dfxp.rds")
z15_end_df = readRDS("data/dfxp/z15_end_df_dfxp.rds")
z15_opc_df = readRDS("data/dfxp/z15_opc_df_dfxp.rds")

d15_oli_df = readRDS("data/dfxp/d15_oli_df_dfxp.rds")
d15_ast_df = readRDS("data/dfxp/d15_ast_df_dfxp.rds")
d15_neu_df = readRDS("data/dfxp/d15_neu_df_dfxp.rds")
d15_mic_df = readRDS("data/dfxp/d15_mic_df_dfxp.rds")
d15_end_df = readRDS("data/dfxp/d15_end_df_dfxp.rds")
d15_opc_df = readRDS("data/dfxp/d15_opc_df_dfxp.rds")

z16_oli_df = readRDS("data/dfxp/z16_oli_df_dfxp.rds")
z16_ast_df = readRDS("data/dfxp/z16_ast_df_dfxp.rds")
z16_neu_df = readRDS("data/dfxp/z16_neu_df_dfxp.rds")
z16_mic_df = readRDS("data/dfxp/z16_mic_df_dfxp.rds")
z16_end_df = readRDS("data/dfxp/z16_end_df_dfxp.rds")

tasic_oli_df = readRDS("data/dfxp/tasic_oli_df_dfxp.rds")
tasic_ast_df = readRDS("data/dfxp/tasic_ast_df_dfxp.rds")
tasic_neu_df = readRDS("data/dfxp/tasic_neu_df_dfxp.rds")
tasic_mic_df = readRDS("data/dfxp/tasic_mic_df_dfxp.rds")
tasic_end_df = readRDS("data/dfxp/tasic_end_df_dfxp.rds")
tasic_opc_df = readRDS("data/dfxp/tasic_opc_df_dfxp.rds")

zeis_oli_df = readRDS("data/dfxp/zeis_oli_df_dfxp.rds")
zeis_ast_df = readRDS("data/dfxp/zeis_ast_df_dfxp.rds")
zeis_neu_df = readRDS("data/dfxp/zeis_neu_df_dfxp.rds")
zeis_mic_df = readRDS("data/dfxp/zeis_mic_df_dfxp.rds")
zeis_end_df = readRDS("data/dfxp/zeis_end_df_dfxp.rds")

sharma_oli_df = readRDS("data/dfxp/sharma_oli_df_dfxp.rds")
sharma_ast_df = readRDS("data/dfxp/sharma_ast_df_dfxp.rds")
sharma_neu_df = readRDS("data/dfxp/sharma_neu_df_dfxp.rds")
sharma_mic_df = readRDS("data/dfxp/sharma_mic_df_dfxp.rds")
sharma_opc_df = readRDS("data/dfxp/sharma_opc_df_dfxp.rds")

#construct volcano plots

make_volcano_facet <- function(data1){

  res_oli = get(paste0(data1, "_oli_df"))
  res_oli$cell = rep("Oligodendrocyte", nrow(res_oli))
  res_mic = get(paste0(data1, "_mic_df"))
  res_mic$cell = rep("Microglia", nrow(res_mic))
  res_neu = get(paste0(data1, "_neu_df"))
  res_neu$cell = rep("Neuron", nrow(res_neu))
  res_ast = get(paste0(data1, "_ast_df"))
  res_ast$cell = rep("Astrocyte", nrow(res_ast))
  data_df = (res_oli, res_mic, res_neu, res_ast)

  if(data1 != "sharma"){
    res_end = get(paste0(data1, "_end_df"))
    res_end$cell = rep("Endothelial", nrow(res_end))
    data_df = rbind(data_df, res_end)
  }

  if(!data1 %in% c("z16", "zeis")){
    res_opc = get(paste0(data1, "_opc_df"))
    res_opc$cell = rep("OPC", nrow(res_opc))
    data_df = rbind(data_df, res_opc)
  }

  data_df$threshold = data_df$logFC > 2 & data_df$adj.P.Val < 0.05

  data1 = data_names_switch(data1)

  volcano = ggplot(data = data_df,
    aes(x = logFC, y = -log10(adj.P.Val), colour = threshold, facet = cell)) +
    geom_point(alpha = 0.95, size = 0.5) +
    facet_wrap( ~ cell, ncol = 2) +
    scale_colour_manual(values = c("#0072B2", "#D55E00"), guide = FALSE) +
    theme_bw() +
    #xlab("Log2 Fold Change") + ylab("-Log10 P-Value")
    ggtitle(data1) + #xlab("") + ylab("") +
    labs(x=NULL, y=NULL) +
    theme(text = element_text(size = 6), strip.background = element_blank()) +  #
    theme(strip.text = element_text(size = 8 , lineheight=0.1)) +
    theme(plot.title = element_text(size = 12))
  # volcano <- ggplotGrob(volcano)
  # volcano$heights[[3]] = unit(0.002, "in")
  #
  return(volcano)

}

z16_volcano = make_volcano_facet("z16")
d15_volcano = make_volcano_facet("d15")
z15_volcano = make_volcano_facet("z15")
sharma_volcano = make_volcano_facet("sharma")
tasic_volcano = make_volcano_facet("tasic")
zeis_volcano = make_volcano_facet("zeis")

grid.arrange(z16_volcano, d15_volcano, z15_volcano,
  sharma_volcano, tasic_volcano, zeis_volcano, ncol = 3)

compare_data_sets <- function(data1, data2, nTerms, orderByCol, orderDecr = TRUE, log_fe = TRUE){

  res_oli = compare_two_signatures(
    get(paste0(data1, "_oli_df")), get(paste0(data2, "_oli_df")), cell = "oli",
    nTerms = nTerms, orderByCol = orderByCol, orderDecr = orderDecr)
  res_mic = compare_two_signatures(
    get(paste0(data1, "_mic_df")), get(paste0(data2, "_mic_df")), cell = "mic",
    nTerms = nTerms, orderByCol = orderByCol, orderDecr = orderDecr)
  res_neu = compare_two_signatures(
    get(paste0(data1, "_neu_df")), get(paste0(data2, "_neu_df")), cell = "neu",
    nTerms = nTerms, orderByCol = orderByCol, orderDecr = orderDecr)
  res_ast = compare_two_signatures(
    get(paste0(data1, "_ast_df")), get(paste0(data2, "_ast_df")), cell = "ast",
    nTerms = nTerms, orderByCol = orderByCol, orderDecr = orderDecr)

  fe_df = rbind(
    res_oli[ , c("fe", "nTerms", "cell")],
    res_neu[ , c("fe", "nTerms", "cell")],
    res_ast[ , c("fe", "nTerms", "cell")],
    res_mic[ , c("fe", "nTerms", "cell")])

  #because sharma doesn't have endothelial cell data
  if(data1 != "sharma" & data2 != "sharma"){
    res_end = compare_two_signatures(
      get(paste0(data1, "_end_df")), get(paste0(data2, "_end_df")), cell = "end",
      nTerms = nTerms, orderByCol = orderByCol, orderDecr = orderDecr)
    fe_df = rbind(fe_df, res_end[ , c("fe", "nTerms", "cell")])
  }

  #because z16 doesn't have OPC cell data
  if((!data1 %in% c("z16", "zeis")) & (!data2 %in% c("z16", "zeis"))){
    res_opc = compare_two_signatures(
      get(paste0(data1, "_opc_df")), get(paste0(data2, "_opc_df")), cell = "opc",
      nTerms = nTerms, orderByCol = orderByCol, orderDecr = orderDecr)
    fe_df = rbind(fe_df, res_opc[ , c("fe", "nTerms", "cell")])
  }

  data1 = data_names_switch(data1)
  data2 = data_names_switch(data2)

  if(log_fe){ fe_df$fe = log(fe_df$fe , 2) }

  print(fe_df)

  fe_layers_plot = ggplot(data = fe_df, aes(x = nTerms, y = fe,
	    group = cell, colour = cell)) + geom_line() + geom_point() +
    theme_bw() +
    # xlab("Number of Top Cell-Type Specific Genes") +
    # ylab("Fold Enrichment") +
    scale_colour_manual("", values = cbPalette, guide = FALSE) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    # guides(fill = guide_legend(title = "Cell")) + #ylim(c(0, 10)) +
    scale_x_log10(breaks = nTerms, labels = nTerms) +
    ggtitle(paste0(data1, " vs. ", data2)) + xlab("") + ylab("") +
    theme(text = element_text(size = 6))

  if(log_fe) { fe_layers_plot = fe_layers_plot + ylim(c(0, 13))}

  fe_df100 = fe_df[fe_df$nTerms == 100, ]

  return(list(plot = fe_layers_plot, fe_df100 = fe_df100))

}

data_names_switch <- function(data_names){

  data_names = gsub("z16", "Zhang (2016)", data_names)
  data_names = gsub("d15", "Darmanis (2015)", data_names)
  data_names = gsub("z15", "Zhang (2015)", data_names)
  data_names = gsub("zeis", "Zeisel (2015)", data_names)
  data_names = gsub("tasic", "Tasic (2016)", data_names)
  data_names = gsub("sharma", "Sharma (2015)", data_names)

  return(data_names)

}

terms_vector = c(10, 20, 50, 100, 200, 500, 1000)

z16_d15 = compare_data_sets(data1 = "z16", data2 = "d15", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
z16_z15 = compare_data_sets(data1 = "z16", data2 = "z15", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
z16_sharma = compare_data_sets(data1 = "z16", data2 = "sharma", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
z16_zeis = compare_data_sets(data1 = "z16", data2 = "zeis", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
z16_tasic = compare_data_sets(data1 = "z16", data2 = "tasic", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)

d15_z15 = compare_data_sets(data1 = "d15", data2 = "z15", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
d15_sharma = compare_data_sets(data1 = "d15", data2 = "sharma", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
d15_zeis = compare_data_sets(data1 = "d15", data2 = "zeis", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
d15_tasic = compare_data_sets(data1 = "d15", data2 = "tasic", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)

z15_sharma = compare_data_sets(data1 = "z15", data2 = "sharma", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
z15_zeis = compare_data_sets(data1 = "z15", data2 = "zeis", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
z15_tasic = compare_data_sets(data1 = "z15", data2 = "tasic", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)

sharma_zeis = compare_data_sets(data1 = "sharma", data2 = "zeis", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)
sharma_tasic = compare_data_sets(data1 = "sharma", data2 = "tasic", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)

zeis_tasic = compare_data_sets(data1 = "zeis", data2 = "tasic", nTerms = terms_vector, orderByCol = "fc_zscore", orderDecr = TRUE)

# > zeis_end_df_match = zeis_end_df[zeis_end_df$genes %in% tasic_end_df$genes, ]
# > tasic_end_df_match = tasic_end_df[tasic_end_df$genes %in% zeis_end_df$genes, ]
# > intersect(head(zeis_end_df, 10)$genes, head(tasic_end_df, 10)$genes )

source("/Users/amckenz/Documents/github/brain_gene_expression/multiplot.R")
grid.arrange(z16_d15[[1]], z16_z15[[1]], z16_sharma[[1]], z16_zeis[[1]],
  z16_tasic[[1]], d15_z15[[1]], d15_sharma[[1]], d15_zeis[[1]], d15_tasic[[1]],
  z15_sharma[[1]], z15_zeis[[1]], z15_tasic[[1]], sharma_zeis[[1]], sharma_tasic[[1]],
  zeis_tasic[[1]], ncol = 3)

# probably just going to ignore the proteomics data for now ...
#z16_sharmap = compare_data_sets(data1 = "z16", data2 = "sharmap", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)

#find the FET p-value at nTerms = 100 (or 200) for all of the pairwise comparisons of data sets, plot them using corrplot for all five major cell types

#identify consensus signatures for all human, mouse, and human + mouse combo

#figure 1: volcano plots (each very small) for each cell type in each data set

#figure 2: pairwise FE of each of the data sets with each other ...

#figure 3: FET p-values for all of the data sets with each other

#supplementary tables: consensus signatures across all of the data sets

data_sets = c("z16", "d15", "z15", "sharma", "zeis", "tasic")
n_data = length(data_sets)
fe_data_oli = matrix(0, nrow = n_data, ncol = n_data)
rownames(fe_data_oli) = colnames(fe_data_oli) = data_sets
fe_data_neu = fe_data_ast = fe_data_mic = fe_data_end = fe_data_opc = fe_data_oli

for(i in 1:n_data){
  for(j in 1:n_data){
    if(i >= j){
      next
    }
    tmp_oli = get(paste0(data_sets[i], "_", data_sets[j]))[["fe_df100"]]
    str(tmp_oli)
    str(tmp_oli[tmp_oli$cell == "oli", "fe"])
    fe_data_oli[i, j] = tmp_oli[tmp_oli$cell == "oli", "fe"]
    tmp_neu = get(paste0(data_sets[i], "_", data_sets[j]))[["fe_df100"]]
    fe_data_neu[i, j] = tmp_neu[tmp_neu$cell == "neu", "fe"]
    tmp_ast = get(paste0(data_sets[i], "_", data_sets[j]))[["fe_df100"]]
    fe_data_ast[i, j] = tmp_ast[tmp_ast$cell == "ast", "fe"]
    tmp_mic = get(paste0(data_sets[i], "_", data_sets[j]))[["fe_df100"]]
    fe_data_mic[i, j] = tmp_mic[tmp_mic$cell == "mic", "fe"]
    if(data_sets[i] == "sharma" | data_sets[j] == "sharma"){
      print("sharma skip")
    } else {
      tmp_end = get(paste0(data_sets[i], "_", data_sets[j]))[["fe_df100"]]
      fe_data_end[i, j] = tmp_end[tmp_end$cell == "end", "fe"]
    }
    if(data_sets[i] %in% c("z16", "zeis")) | data_sets[j] %in% c("z16", "zeis"))){
      print("z16 skip")
    } else {
      tmp_opc = get(paste0(data_sets[i], "_", data_sets[j]))[["fe_df100"]]
      str(tmp_opc)
      fe_data_opc[i, j] = tmp_opc[tmp_opc$cell == "opc", "fe"]
    }
  }
}

# m[lower.tri(m)] = t(m)[lower.tri(m)]

fe_data_oli[lower.tri(fe_data_oli)] = t(fe_data_oli)[lower.tri(fe_data_oli)]
corrplot.mixed(fe_data_oli, lower = "number", upper = "circle", tl.pos = "lt", addgrid.col = "black", is.corr = FALSE)

fe_data_neu[lower.tri(fe_data_neu)] = t(fe_data_neu)[lower.tri(fe_data_neu)]
corrplot.mixed(fe_data_neu, lower = "number", upper = "circle", tl.pos = "lt", addgrid.col = "black", is.corr = FALSE)

fe_data_mic[lower.tri(fe_data_mic)] = t(fe_data_mic)[lower.tri(fe_data_mic)]
corrplot.mixed(fe_data_mic, lower = "number", upper = "circle", tl.pos = "lt", addgrid.col = "black", is.corr = FALSE)

fe_data_ast[lower.tri(fe_data_ast)] = t(fe_data_ast)[lower.tri(fe_data_ast)]
corrplot.mixed(fe_data_ast, lower = "number", upper = "circle", tl.pos = "lt", addgrid.col = "black", is.corr = FALSE)

fe_data_end = fe_data_end[!rownames(fe_data_end) %in% "sharma", !colnames(fe_data_end) %in% "sharma"]
fe_data_end[lower.tri(fe_data_end)] = t(fe_data_end)[lower.tri(fe_data_end)]
corrplot.mixed(fe_data_end,  lower = "number", upper = "circle", tl.pos = "lt", addgrid.col = "black", is.corr = FALSE)

fe_data_opc = fe_data_opc[!rownames(fe_data_opc) %in% c("z16", "zeis"), !colnames(fe_data_opc) %in% c("z16", "zeis")]
fe_data_opc[lower.tri(fe_data_opc)] = t(fe_data_opc)[lower.tri(fe_data_opc)]
corrplot.mixed(fe_data_opc, lower = "number", upper = "circle", tl.pos = "lt", addgrid.col = "black", is.corr = FALSE)

get_forest_plot <- function(gene, cell, auto_min = FALSE){

  d15_res = get(paste0("d15_", cell, "_df"))
  d15_res$study = rep("Darmanis (2015)", nrow(d15_res))
  d15_res$type = rep("Human", nrow(d15_res))
  z15_res = get(paste0("z15_", cell, "_df"))
  z15_res$study = rep("Zhang (2015)", nrow(z15_res))
  z15_res$type = rep("Mouse", nrow(z15_res))
  tasic_res = get(paste0("tasic_", cell, "_df"))
  tasic_res$study = rep("Tasic (2016)", nrow(tasic_res))
  tasic_res$type = rep("Mouse", nrow(tasic_res))
  cell_df = rbind(d15_res, z15_res, tasic_res)

  if(cell != "opc"){
    zeis_res = get(paste0("zeis_", cell, "_df"))
    zeis_res$study = rep("Zeisel (2015)", nrow(zeis_res))
    zeis_res$type = rep("Mouse", nrow(zeis_res))
    z16_res = get(paste0("z16_", cell, "_df"))
    z16_res$study = rep("Zhang (2016)", nrow(z16_res))
    z16_res$type = rep("Human", nrow(z16_res))
    cell_df = rbind(cell_df, zeis_res, z16_res)
  }

  if(cell != "end"){
    sharma_res = get(paste0("sharma_", cell, "_df"))
    sharma_res$study = rep("Sharma (2015)", nrow(sharma_res))
    sharma_res$type = rep("Mouse", nrow(sharma_res))
    cell_df = rbind(cell_df, sharma_res)
  }

  cell_df = cell_df[cell_df$genes %in% gene, ]

  print(cell_df)

  forest_plot = ggplot(cell_df, aes(x = logFC, y = study, xmin = CI.L, xmax = CI.R, color = type)) +
    geom_vline(xintercept = 0.0, linetype = 2, alpha = 0.75) +
    geom_errorbarh(alpha = 0.5, height = 0) +
    geom_point(shape = 15) +
    theme_bw() + ggtitle(paste0(gene, " in ", toupper(cell))) + ylab(NULL) + xlab("log2 Fold-Change") +
    theme(legend.position="none") +
    theme(text = element_text(size = 8)) +  #
    theme(plot.title = element_text(size = 12))

  if(auto_min) { forest_plot = forest_plot + xlim(min(0, cell_df$CI.L), max(cell_df$CI.R)) }
  if(!auto_min) { forest_plot = forest_plot + xlim(-1.75, 17.5) }


  return(forest_plot)

}

plp1 = get_forest_plot("PLP1", "oli")
snap25 = get_forest_plot("SNAP25", "neu")
ctss = get_forest_plot("CTSS", "mic")
itm2a = get_forest_plot("ITM2A", "end")
gja1 = get_forest_plot("GJA1", "ast")
pdgfra = get_forest_plot("PDGFRA", "opc")
grid.arrange(plp1, snap25, gja1, ctss, itm2a, pdgfra, ncol = 2)

find_deg_number <- function(data1){

  res_oli = get(paste0(data1, "_oli_df"))
  res_oli$cell = rep("Oligodendrocyte", nrow(res_oli))
  n_deg_oli = sum(res_oli$logFC > 2 & res_oli$adj.P.Val < 0.05)
  res_mic = get(paste0(data1, "_mic_df"))
  res_mic$cell = rep("Microglia", nrow(res_mic))
  n_deg_mic = sum(res_mic$logFC > 2 & res_mic$adj.P.Val < 0.05)
  res_neu = get(paste0(data1, "_neu_df"))
  res_neu$cell = rep("Neuron", nrow(res_neu))
  n_deg_neu = sum(res_neu$logFC > 2 & res_neu$adj.P.Val < 0.05)
  res_ast = get(paste0(data1, "_ast_df"))
  res_ast$cell = rep("Astrocyte", nrow(res_ast))
  n_deg_ast = sum(res_ast$logFC > 2 & res_ast$adj.P.Val < 0.05)
  n_deg_vec = rbind(n_deg_oli, n_deg_mic, n_deg_neu, n_deg_ast)

  if(data1 != "sharma"){
    res_end = get(paste0(data1, "_end_df"))
    res_end$cell = rep("Endothelial", nrow(res_end))
    n_deg_end = sum(res_end$logFC > 2 & res_end$adj.P.Val < 0.05)
    n_deg_vec = rbind(n_deg_vec, n_deg_end)
  }

  if(!data1 %in% c("z16", "zeis")){
    res_opc = get(paste0(data1, "_opc_df"))
    res_opc$cell = rep("OPC", nrow(res_opc))
    n_deg_opc = sum(res_opc$logFC > 2 & res_opc$adj.P.Val < 0.05)
    n_deg_vec = rbind(n_deg_vec, n_deg_opc)
  }

  return(n_deg_vec)

}

std_error <- function(x) sd(x)/sqrt(length(x))
full_deg_res = as.numeric(find_deg_number("z16"), find_deg_number("d15"), find_deg_number("z15"),
  find_deg_number("zeis"), find_deg_number("tasic"), find_deg_number("sharma"))
mean(full_deg_res)
std_error(full_deg_res)

plot_set <- function(cell, nGenes){

  d15_res = get(paste0("d15_", cell, "_df"))
  z15_res = get(paste0("z15_", cell, "_df"))
  tasic_res = get(paste0("tasic_", cell, "_df"))
  gene_sets = list(d15_res$genes, z15_res$genes, tasic_res$genes)

  if(cell != "opc"){
    zeis_res = get(paste0("zeis_", cell, "_df"))
    z16_res = get(paste0("z16_", cell, "_df"))
    gene_sets$a = z16_res$genes
    gene_sets$b = zeis_res$genes
  }

  if(cell != "end"){
    sharma_res = get(paste0("sharma_", cell, "_df"))
    gene_sets$c = sharma_res$genes
  }

  universe = Reduce(intersect, gene_sets)
  universe = unique(universe)

  d15_res = d15_res[d15_res$genes %in% universe, ]
  z15_res = z15_res[z15_res$genes %in% universe, ]
  tasic_res = tasic_res[tasic_res$genes %in% universe, ]

  gene_sets = list(
    head(d15_res$genes, nGenes),
    head(z15_res$genes, nGenes),
    head(tasic_res$genes, nGenes))

  data_sets = c("d15", "z15", "tasic")

  if(cell != "opc"){
    zeis_res = zeis_res[zeis_res$genes %in% universe, ]
    z16_res = z16_res[z16_res$genes %in% universe, ]
    gene_sets$a = head(z16_res$genes, nGenes)
    gene_sets$b = head(zeis_res$genes, nGenes)
    data_sets = c(data_sets, "z16", "zeis")
  }

  if(cell != "end"){
    sharma_res = sharma_res[sharma_res$genes %in% universe, ]
    gene_sets$c = head(sharma_res$genes, nGenes)
    data_sets = c(data_sets, "sharma")
  }

  data_names = data_names_switch(data_sets)
  str(data_names)
  names(gene_sets) = data_names
  str(gene_sets)

  cell_res = supertest(gene_sets, n = length(universe))
  cell_res$P.value = p.adjust(cell_res$P.value, method = "BH")

  cellres_df = summary(cell_res)$Table
  cellres_df = cellres_df[!cellres_df$Degree == 1,]
  cellres_df = cellres_df[ , !colnames(cellres_df) %in% "Elements"]
  # cellres_df$adj_pval = p.adjust(cellres_df$P.value, method = "BH")

  # print(cellres_df)
  print(cellres_df[cellres_df$Intersections == "Darmanis (2015) & Zhang (2016)", ])
  print(cellres_df[cellres_df$Intersections == "Zhang (2015) & Tasic (2016) & Zeisel (2015) & Sharma (2015)", ])
  print(cellres_df[cellres_df$Intersections == "Darmanis (2015) & Zhang (2015) & Tasic (2016)", ])

  plot(cell_res, sort.by = 'size', degree = 2:length(gene_sets),
    keep.empty.intersections = FALSE, cex = 1, legend.pos=c(2,2)) #circular


}

plot_set("oli", 200)
plot_set("neu", 200)
plot_set("mic", 200)
plot_set("ast", 200)
plot_set("end", 200)
plot_set("opc", 200)


#this function ranks the genes across multiple data sets
rank_genes_logfc <- function(cell, subset = NULL, mouse_vs_human = FALSE, orderDown = TRUE){

  drops = NULL
  if(cell == "opc"){ drops = c(drops, "z16", "zeis") }
  if(cell == "end"){ drops = c(drops, "sharma") }
  if(!is.null(subset)){
    if(subset == "mouse"){ drops = c(drops, "z16", "d15") }
    if(subset == "human"){ drops = c(drops, "z15", "sharma", "zeis", "tasic") }
  }
  if(length(drops) == 5){ stop("No consensus is available for this cell type in this subset because there is only one data set that fit the criteria.") }

  data_sets = c("z16", "d15", "z15", "sharma", "zeis", "tasic")
  if(!is.null(drops)){
    data_keep = data_sets[!data_sets %in% drops]
  } else {
    data_keep = data_sets
  }

  n_data_sets = length(data_keep)

  #merge the data sets into one, keep all gene symbols
  for(i in 1:n_data_sets){
    data_keep[i]
    tmp = get(paste0(data_keep[i], "_", cell, "_df"))
    tmp$genes = toupper(tmp$genes)
    for(j in 2:(ncol(tmp) - 1)){
      colnames(tmp)[j] = paste0(colnames(tmp)[j], "_", data_keep[i])
    }
    if(i == 1){
      merged = tmp
    } else {
      merged = merge(merged, tmp, by = "genes", all = TRUE) #, all = TRUE
    }
  }

  #filter to remove genes that are present in too few of the data sets
  log_fc_only_merged = merged[ , grepl("fc_zscore", colnames(merged))]
  merged$n_data_na = apply(log_fc_only_merged, 1, function(x) sum(is.na(x)))
  merged = merged[!merged$n_data_na >= 0.5, ]
  #if merging mice and human, require at least one of each
  if(is.null(subset)){
    log_fc_only_merged = merged[ , grepl("fc_zscore", colnames(merged))]
    merged$non_data_na_human = apply(log_fc_only_merged, 1, function(x) any(!is.na(x[data_sets %in% c("z16", "d15")])))
    merged$non_data_na_mice = apply(log_fc_only_merged, 1, function(x) any(!is.na(x[data_sets %in% c("z15", "sharma", "zeis", "tasic")])))
    merged = merged[merged$non_data_na_human, ]
    merged = merged[merged$non_data_na_mice, ]
  }

  #calculate the grand mean
  log_fc_merged = merged[ , grepl("fc_zscore", colnames(merged))]
  log_fc_merged$grand_mean = apply(log_fc_merged, 1, median, na.rm = TRUE)
  log_fc_merged$gene = merged$genes
  log_fc_merged = log_fc_merged[order(log_fc_merged$grand_mean, decreasing = TRUE), ]

  #mouse vs human difference
  if(mouse_vs_human){
    log_fc_merged = merged[ , grepl("fc_zscore", colnames(merged))]
    fc_res = vector(mode = "numeric", length = nrow(log_fc_merged))
    t_res = vector(mode = "numeric", length = nrow(log_fc_merged))
    p_vals = vector(mode = "numeric", length = nrow(log_fc_merged))
    for(i in 1:nrow(log_fc_merged)){
      str(as.numeric(log_fc_merged[i, data_sets %in% c("z16", "d15")]))
      str(as.numeric(log_fc_merged[i, data_sets %in% c("z15", "sharma", "zeis", "tasic")]))
      fc_res[i] = mean(as.numeric(log_fc_merged[i, data_sets %in% c("z16", "d15")]), na.rm = TRUE) -
        mean(as.numeric(log_fc_merged[i, data_sets %in% c("z15", "sharma", "zeis", "tasic")]), na.rm = TRUE)
      t_res[i] = t.test(as.numeric(log_fc_merged[i, data_sets %in% c("z16", "d15")]),
        as.numeric(log_fc_merged[i, data_sets %in% c("z15", "sharma", "zeis", "tasic")]), na.rm = TRUE)$statistic
      p_vals[i] = t.test(as.numeric(log_fc_merged[i, data_sets %in% c("z16", "d15")]),
        as.numeric(log_fc_merged[i, data_sets %in% c("z15", "sharma", "zeis", "tasic")]), na.rm = TRUE)$p.value
    }
    log_fc_merged$mouse_vs_human_fc_res = fc_res
    log_fc_merged$mouse_vs_human_t_res = t_res
    log_fc_merged$mouse_vs_human_p_vals = p_vals
    log_fc_merged$mouse_vs_human_p_val_adj = qvalue(p_vals)$qvalue
    log_fc_merged$gene = merged$genes
    # log_fc_merged$mouse_vs_human_res =
    #   mean(as.numeric(log_fc_merged[, data_sets %in% c("z16", "d15")])) -
    #   mean(as.numeric(log_fc_merged[, data_sets %in% c("z15", "sharma", "zeis", "tasic")]))
    log_fc_merged = log_fc_merged[order(log_fc_merged$mouse_vs_human_p_vals, decreasing = orderDown), ]
    return(log_fc_merged)
  } else {
    return(log_fc_merged)
  }

}

neu_top = rank_genes_logfc("neu")
mic_top = rank_genes_logfc("mic")
ast_top = rank_genes_logfc("ast")
end_top = rank_genes_logfc("end")
opc_top = rank_genes_logfc("opc")
oli_top = rank_genes_logfc("oli")

write.delim(neu_top, "neu_top_all.tsv")
write.delim(mic_top, "mic_top_all.tsv")
write.delim(ast_top, "ast_top_all.tsv")
write.delim(end_top, "end_top_all.tsv")
write.delim(opc_top, "opc_top_all.tsv")
write.delim(oli_top, "oli_top_all.tsv")

oli_top_human = rank_genes_logfc("oli", subset = "human")
neu_top_human = rank_genes_logfc("neu", subset = "human")
ast_top_human = rank_genes_logfc("ast", subset = "human")
mic_top_human = rank_genes_logfc("mic", subset = "human")
end_top_human = rank_genes_logfc("end", subset = "human")
opc_top_human = rank_genes_logfc("opc", subset = "human")

write.delim(oli_top_human, "oli_top_human.tsv")
write.delim(neu_top_human, "neu_top_human.tsv")
write.delim(ast_top_human, "ast_top_human.tsv")
write.delim(mic_top_human, "mic_top_human.tsv")
write.delim(end_top_human, "end_top_human.tsv")
write.delim(opc_top_human, "opc_top_human.tsv")
