
library(HGNChelper)
library(ggplot2)
library(DGCA)
library(bayesbio)

#reads in

source("/Users/amckenz/Dropbox/zhang/general_code/get_fold_enrichment.R")

#a function that takes two of the df's, orders and identifies top N with one of the cutoff types,
#then find jaccard index of the plot, and spit back the results.
#part of this function needs to clean the gene symbols with HGNChelper
#and then purge/merge the df's so that only the shared genes are included.
compare_two_signatures <- function(df1, df2, nTerms, orderByCol, orderDecr, cell){

  df1$Row.names = switchGenesToHGCN(df1$Row.names)
  df2$Row.names = switchGenesToHGCN(df2$Row.names)
  universe = intersect(df1$Row.names, df2$Row.names)
  df1 = df1[df1$Row.names %in% universe, ]
  df2 = df2[df2$Row.names %in% universe, ]
  str(df1)
  str(df2)
  df2 = df2[!duplicated(df2$Row.names), ]
  str(df2)
  str(universe)

  df1 = df1[order(df1[ , orderByCol], decreasing = orderDecr), ]
  df2 = df2[order(df2[ , orderByCol], decreasing = orderDecr), ]

  # jac_res = vector()
  # for(i in 1:length(nTerms)){
  #   top_df1 = head(df1, nTerms[i])$Row.names
  #   top_df2 = head(df2, nTerms[i])$Row.names
  #   jac_res[i]= jaccardSets(top_df1, top_df2)
  # }
  # return(jac_res)

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

terms_vector = c(10, 20, 50, 100, 200, 500, 1000)

# res = compare_two_signatures(z16_mic_df, d15_mic_df, nTerms = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 15000), orderByCol = "mean_fc", orderDecr = FALSE)
# res = compare_two_signatures(z16_mic_df, d15_mic_df, nTerms = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000), orderByCol = "logFC", orderDecr = FALSE)
res_mic = compare_two_signatures(z16_mic_df, d15_mic_df, nTerms = terms_vector, orderByCol = "P.Value", orderDecr = FALSE, cell = "mic")

# res = compare_two_signatures(z16_oli_df, d15_oli_df, nTerms = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000), orderByCol = "mean_fc", orderDecr = FALSE)
# res = compare_two_signatures(z16_oli_df, d15_oli_df, nTerms = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000), orderByCol = "logFC", orderDecr = FALSE)
res_oli = compare_two_signatures(z16_oli_df, d15_oli_df, nTerms = terms_vector, orderByCol = "P.Value", orderDecr = FALSE, cell = "oli")

res_neu = compare_two_signatures(z16_neu_df, d15_neu_df, nTerms = terms_vector, orderByCol = "P.Value", orderDecr = FALSE, cell = "neu")
res_ast = compare_two_signatures(z16_ast_df, d15_ast_df, nTerms = terms_vector, orderByCol = "P.Value", orderDecr = FALSE, cell = "ast")
res_end = compare_two_signatures(z16_end_df, d15_end_df, nTerms = terms_vector, orderByCol = "P.Value", orderDecr = FALSE, cell = "end")



# intersect(head(z16_mic_df, 50)$Row.names, head(d15_mic_df, 50)$Row.names)

#a function that combines the jaccard indices from each cell type into a plot
#with like 10, 20, 50, 100, 200, 500, 1000 top genes from each ...
# plot_jaccard_log_terms <- function(){
# }

fe_df = rbind(
  res_mic[ , c("fe", "nTerms", "cell")],
  res_oli[ , c("fe", "nTerms", "cell")],
  res_neu[ , c("fe", "nTerms", "cell")],
  res_ast[ , c("fe", "nTerms", "cell")],
  res_end[ , c("fe", "nTerms", "cell")])

cbPalette = c("#CC79A7", "#E69F00", "#56B4E9", "#000000", "#009E73", "#F0E442", "#0072B2",
 "red", "lightgreen", "#999999", "#990000")

fe_layers_plot = ggplot(data = fe_df, aes(x = factor(nTerms), y = log(fe, 2),
	    group = cell, colour = cell)) + geom_line() + geom_point() +
    theme_bw() + xlab("Number of Top Cell-Type Specific Genes") +
    ylab("Fold Enrichment") + scale_colour_manual("", values = cbPalette) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    guides(fill = guide_legend(title = "Cell"))
