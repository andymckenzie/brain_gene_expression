library(qvalue)

shrink_pval_and_avg_expr <- function(degs){

  locfdr = qvalue::qvalue(degs$P.Value)$lfdr

  #formula from 2.2.3 http://biostats.bepress.com/cobra/art60/;
  #also see http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-63
  pval_shrink_logFC = (1 - locfdr) * degs$logFC

  #formula from eq. 1 and 2 http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2464587/
  shrunken_logFC = pval_shrink_logFC *
    ((degs$expr_average - min(degs$expr_average)) /
    (max(degs$expr_average) - min(degs$expr_average)))

  #convert to z-scores
  fc_zscore = scale(shrunken_logFC, center = TRUE, scale = TRUE)

  return(fc_zscore)

}
