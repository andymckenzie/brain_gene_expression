
library(HGNChelper)
library(ggplot2)
library(DGCA)
library(bayesbio)
library(corrplot)

cbPalette = c("#CC79A7", "#E69F00", "#56B4E9", "#000000", "#009E73", "#F0E442", "#0072B2",
 "red", "lightgreen", "#999999", "#990000")

setwd("/Users/amckenz/Documents/github/brain_gene_expression/")

source("/Users/amckenz/Dropbox/zhang/general_code/get_fold_enrichment.R")

#a function that takes two of the df's, orders and identifies top N with one of the cutoff types,
#then find jaccard index of the plot, and spit back the results.
#part of this function needs to clean the gene symbols with HGNChelper
#and then purge/merge the df's so that only the shared genes are included.
compare_two_signatures <- function(df1, df2, nTerms, orderByCol, orderDecr, cell, HGNC_clean = FALSE, toupper = TRUE){

  if(HGNC_clean){
    df1$Row.names = switchGenesToHGCN(df1$Row.names)
    df2$Row.names = switchGenesToHGCN(df2$Row.names)
  }
  if(toupper){
    df1$Row.names = toupper(df1$Row.names)
    df2$Row.names = toupper(df2$Row.names)
  }
  universe = intersect(df1$Row.names, df2$Row.names)
  df1 = df1[df1$Row.names %in% universe, ]
  df2 = df2[df2$Row.names %in% universe, ]
  df2 = df2[!duplicated(df2$Row.names), ]

  df1 = df1[order(df1[ , orderByCol], decreasing = orderDecr), ]
  df2 = df2[order(df2[ , orderByCol], decreasing = orderDecr), ]

  fe = vector()
  pvals = vector()

  for(i in 1:length(nTerms)){
    top_n_df1 = unique(head(df1, nTerms[i])$Row.names)
    top_n_df2 = unique(head(df2, nTerms[i])$Row.names)
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
zeis_opc_df = readRDS("data/dfxp/zeis_opc_df_dfxp.rds")

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
  str(res_oli)
  str(res_mic)
  str(res_neu)
  str(res_ast)
  data_df = rbind(res_oli, res_mic, res_neu, res_ast)

  if(data1 != "sharma"){
    res_end = get(paste0(data1, "_end_df"))
    res_end$cell = rep("Endothelial", nrow(res_end))
    data_df = rbind(data_df, res_end)
  }

  if(data1 != "z16"){
    res_opc = get(paste0(data1, "_opc_df"))
    res_opc$cell = rep("OPC", nrow(res_opc))
    data_df = rbind(data_df, res_opc)
  }

  data_df$threshold = data_df$logFC > 2 & data_df$P.Value < 0.01

  volcano = ggplot(data = data_df,
    aes(x = logFC, y = -log10(adj.P.Val), colour = threshold, facet = cell)) +
    geom_point(alpha = 0.95, size = 0.5) +
    facet_wrap( ~ cell, ncol = 2) +
    scale_colour_manual(values = c("#0072B2", "#D55E00"), guide = FALSE) +
    xlab("Log2 Fold Change") + ylab("-Log10 P-Value") + theme_bw()

  return(volcano)

}

make_volcano_facet("zm16")






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
  if(data1 != "z16" & data2 != "z16"){
    res_opc = compare_two_signatures(
      get(paste0(data1, "_opc_df")), get(paste0(data2, "_opc_df")), cell = "opc",
      nTerms = nTerms, orderByCol = orderByCol, orderDecr = orderDecr)
    fe_df = rbind(fe_df, res_opc[ , c("fe", "nTerms", "cell")])
  }


  if(log_fe){ fe_df$fe = log(fe_df$fe , 2) }

  fe_layers_plot = ggplot(data = fe_df, aes(x = nTerms, y = fe,
	    group = cell, colour = cell)) + geom_line() + geom_point() +
    theme_bw() + xlab("Number of Top Cell-Type Specific Genes") +
    ylab("Fold Enrichment") + scale_colour_manual("", values = cbPalette) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    guides(fill = guide_legend(title = "Cell")) + ylim(c(0, 10)) +
    scale_x_log10(breaks = nTerms, labels = nTerms)

  fe_df100 = fe_df[fe_df$nTerms == 100, ]

  return(list(plot = fe_layers_plot, fe_df100 = fe_df100))

}

terms_vector = c(10, 20, 50, 100, 200, 500, 1000)

z16_d15 = compare_data_sets(data1 = "z16", data2 = "d15", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
z16_z15 = compare_data_sets(data1 = "z16", data2 = "z15", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
z16_sharma = compare_data_sets(data1 = "z16", data2 = "sharma", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
z16_zeis = compare_data_sets(data1 = "z16", data2 = "zeis", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
z16_tasic = compare_data_sets(data1 = "z16", data2 = "tasic", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)

d15_z15 = compare_data_sets(data1 = "d15", data2 = "z15", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
d15_sharma = compare_data_sets(data1 = "d15", data2 = "sharma", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
d15_zeis = compare_data_sets(data1 = "d15", data2 = "zeis", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
d15_tasic = compare_data_sets(data1 = "d15", data2 = "tasic", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)

z15_sharma = compare_data_sets(data1 = "z15", data2 = "sharma", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
z15_zeis = compare_data_sets(data1 = "z15", data2 = "zeis", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
z15_tasic = compare_data_sets(data1 = "z15", data2 = "tasic", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)

sharma_zeis = compare_data_sets(data1 = "sharma", data2 = "zeis", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)
sharma_tasic = compare_data_sets(data1 = "sharma", data2 = "tasic", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)

zeis_tasic = compare_data_sets(data1 = "zeis", data2 = "tasic", nTerms = terms_vector, orderByCol = "logFC", orderDecr = TRUE)

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
    if(data_sets[i] == "z16" | data_sets[j] == "z16"){
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

fe_data_opc = fe_data_opc[!rownames(fe_data_opc) %in% "z16", !colnames(fe_data_opc) %in% "z16"]
fe_data_opc[lower.tri(fe_data_opc)] = t(fe_data_opc)[lower.tri(fe_data_opc)]
corrplot.mixed(fe_data_opc, lower = "number", upper = "circle", tl.pos = "lt", addgrid.col = "black", is.corr = FALSE)

plot_meta_analysis <- function(gene, cell, subset = NULL){

  drops = NULL

  if(cell == "opc"){
    drops = c(drops, "z16")
  }
  if(cell == "end"){
    drops = c(drops, "sharma")
  }



}
