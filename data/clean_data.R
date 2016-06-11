library(readxl)

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

saveRDS(sharma_num, "sharma_rpkm.rds")

#read in the sharma proteomics data
sharma_prot = read_excel("data/sharma_proteomics_nn.4160-S7.xlsx", skip = 1)
sharma_prot_num = sharma_prot[, !colnames(sharma_prot) %in%
  c("Gene names", "Protein names", "PEP", "Mol. weight [kDa]",
  "Sequence coverage [%]", "Protein IDs", "Majority protein IDs")]
sharma_prot$`Gene names`[is.na(sharma_prot$`Gene names`)] = ""
rownames(sharma_prot_num) = make.unique(sharma_prot$`Gene names`)
#proteomics data is also approximately normalized... colMeans(as.matrix(sharma_num))

saveRDS(sharma_prot_num, "sharma_lfq.rds")
