library(limma)

dfxp_function <- function(contrast_matrix, cell_type, list_other_cells){

	expr_average = rowMeans(barres_merged_log[ , which(sample_types %in% cell_type)])
	cutoff = as.numeric(quantile(expr_average, cutoff_thresh, names = FALSE))
	barres_trimmed = barres_merged_log[(expr_average > cutoff), ]
	#print(nrow(barres_trimmed))

	design = model.matrix(~ fl + 0, barres_trimmed)
	colnames(design) = levels(fl)
	fit = lmFit(barres_trimmed, design)

	fit2 = contrasts.fit(fit, contrast_matrix)
	fit2 = eBayes(fit2, 0.01, trend = TRUE)
	tT = topTable(fit2, adjust="none", sort.by="B", number = nrow(fit2))

	#print(grep("Prkd2", row.names(tT)))
	tT = tT[tT$P.Value < pval_thresh, ]

	#compute fold-changes for CTOI vs others
	for(i in 1:length(list_other_cells)){
		fc = rowMeans(barres_trimmed[ , which(sample_types %in% cell_type)]) /
			rowMeans(barres_trimmed[ , which(sample_types %in% list_other_cells[[i]])])
		if(!exists("fc_df")){
			fc_df = data.frame(fc)
		} else {
			fc_df = cbind(fc_df, fc)
		}
	}
	fc_min = apply(fc_df, 1, min)
	fc_min_names = data.frame(row.names(barres_trimmed), fc_min)
	names(fc_min_names) = c("barres_names", "min_FC")

	#since already sorted, need to merge via rownames

	# fc = rowMeans(barres_trimmed[ , which(sample_types %in% cell_type)]) /
	# 	rowMeans(barres_trimmed[ , which(sample_types %nin% cell_type & sample_types %nin% "WC")])

	# print(grep("Prkd2", row.names(tT)))
	# print(tT[grep("Prkd2", row.names(tT)), ])

	toptable = merge(tT, fc_min_names, by.x = "row.names", by.y = "barres_names")

	toptable = toptable[toptable$min_FC > fc_thres, ]

	toptable = toptable[order(toptable$P.Value), ]

	return(toptable)

}
