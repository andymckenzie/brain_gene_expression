library(shiny)
library(ggplot2)

cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

shinyServer(function(input, output) {

  output$cellPlot <- renderPlot(height = 600, units = "px", width = 600, {

  std_error <- function(x) sd(x)/sqrt(length(x))

  if(input$data_type == "Sharma_RNA"){

    get_sharma_rpkm_plot <- function(gene){
      if(!gene %in% rownames(sharma_rpkm)){
        gene = cap_first(gene)
        if(!gene %in% rownames(sharma_rpkm)){
          stop("Even when uncapitalized, that gene symbol is not in the Sharma et al. RPKM data set.")
        }
      }
    	mean = vector()
    	se = vector()
    	for(i in 0:6){
        if(input$log){
      		mean[i+1] = mean(log(as.numeric(sharma_rpkm[gene, (i*3+1):((i*3)+3)]), 2))
      		se[i+1] = std_error(log(as.numeric(sharma_rpkm[gene, (i*3+1):((i*3)+3)]), 2))
        } else {
          mean[i+1] = mean(as.numeric(sharma_rpkm[gene, (i*3+1):((i*3)+3)]))
          se[i+1] = std_error(as.numeric(sharma_rpkm[gene, (i*3+1):((i*3)+3)]))
        }
    	}
    	types = c("Oligodendrocytes Div1", "Oligodendrocytes Div2.5", "Oligodendrocytes Div4",
        "Microglia", "Astrocyte", "Neurons Div0.5", "Neurons Div10")
    	tmp = data.frame(mean = mean, se = se, types = types)
      tmp$types = factor(tmp$types, levels = unique(tmp$types))
    	limits = aes(ymax = mean + se, ymin = mean - se)
    	tmp_plot = ggplot(data = tmp, aes(x = types, y = mean)) + geom_point() +
    		geom_errorbar(limits, width = 0.2) + theme_bw() +
    		ylab("RNA Expression") + xlab("") +
    		theme(axis.text.x = element_text(angle = 60, hjust = 1),
          text = element_text(size = 22))
      print(tmp_plot)
    }

    get_sharma_rpkm_plot(input$gene)

  }

  if(input$data_type == "Sharma_Proteomics"){

    get_sharma_lfq_plot <- function(gene){
      if(!gene %in% rownames(sharma_lfq)){
        gene = cap_first(gene)
        if(!gene %in% rownames(sharma_lfq)){
          stop("Even when uncapitalized, that gene symbol is not in the Sharma et al. LFQ data set.")
        }
      }
      mean = vector()
      se = vector()
      for(i in 0:8){
        if(input$log){
          mean[i+1] = mean(log(as.numeric(sharma_lfq[gene, (i*3+1):((i*3)+3)]), 2), na.rm = TRUE)
          se[i+1] = std_error(log(as.numeric(sharma_lfq[gene, (i*3+1):((i*3)+3)]), 2))
        } else {
          mean[i+1] = mean(as.numeric(sharma_lfq[gene, (i*3+1):((i*3)+3)]))
          se[i+1] = std_error(as.numeric(sharma_lfq[gene, (i*3+1):((i*3)+3)]))
        }
      }
      types = c("Adult Microglia", "Young Microglia", "Astrocytes",
        "Neurons Div5", "Neurons Div10", "Neurons Div15",
        "Oligodendrocytes Div1", "Oligodendrocytes Div2.5", "Oligodendrocytes Div4")
      tmp = data.frame(mean = mean, se = se, types = types)
      tmp$types = factor(tmp$types, levels = unique(tmp$types))
      limits = aes(ymax = mean + se, ymin = mean - se)
      tmp_plot = ggplot(data = tmp, aes(x = types, y = mean)) + geom_point() +
        geom_errorbar(limits, width = 0.2) + theme_bw() +
        ylab("Protein Expression") + xlab("") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
          text = element_text(size = 22))
      print(tmp_plot)
    }

    get_sharma_lfq_plot(input$gene)

  }

  if(input$data_type == "Darmanis_RNA"){

    get_darmanis_plot <- function(gene){

      if(!gene %in% rownames(darmanis)){
        gene = toupper(gene)
        if(!gene %in% rownames(darmanis)){
          stop("Even when capitalized, that gene symbol is not in the Darmanis et al. data set.")
        }
      }
      gene_data = darmanis[gene, ]
      if(input$log){
        types = colnames(darmanis)[grepl("log", colnames(darmanis))]
        tmp = gene_data[ , types]
        mean = tmp[ , grepl("mean", colnames(tmp))]
        se = tmp[ , grepl("se", colnames(tmp))]
      } else {
        types = colnames(darmanis)[!grepl("log", colnames(darmanis))]
        tmp = gene_data[ , types]
        mean = tmp[ , grepl("mean", colnames(tmp))]
        se = tmp[ , grepl("se", colnames(tmp))]
      }
      types_clean = types[grepl("mean", types)]
      types_clean = sapply(strsplit(types_clean, "_", fixed = TRUE), "[", 1)
      tmp_df = data.frame(mean = as.numeric(mean), se = as.numeric(se), types = types_clean)
      tmp_df$types = factor(tmp_df$types, levels = unique(tmp_df$types))
      limits = aes(ymax = mean + se, ymin = mean - se)
      tmp_plot = ggplot(data = tmp_df, aes(x = types, y = mean)) + geom_point() +
        geom_errorbar(limits, width = 0.2) + theme_bw() +
        ylab("TMM-Normalized RNA CPM") + xlab("") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
          text = element_text(size = 22))
      print(tmp_plot)

    }

    get_darmanis_plot(input$gene)

  }

  if(input$data_type == "Marques_RNA"){

    get_marques_plot <- function(gene){

      if(!gene %in% rownames(marques)){
        gene = cap_first(gene)
        if(!gene %in% rownames(marques)){
          stop("Even when uncapitalized, that gene symbol is not in the Marques et al. data set.")
        }
      }
      gene_data = marques[gene, ]
      if(input$log){
        types = colnames(marques)[grepl("log", colnames(marques))]
        tmp = gene_data[ , types]
        mean = tmp[ , grepl("mean", colnames(tmp))]
        se = tmp[ , grepl("se", colnames(tmp))]
      } else {
        types = colnames(marques)[!grepl("log", colnames(marques))]
        tmp = gene_data[ , types]
        mean = tmp[ , grepl("mean", colnames(tmp))]
        se = tmp[ , grepl("se", colnames(tmp))]
      }
      types_clean = types[grepl("mean", types)]
      types_clean = sapply(strsplit(types_clean, "_", fixed = TRUE), "[", 1)
      tmp_df = data.frame(mean = as.numeric(mean), se = as.numeric(se), types = types_clean)
      ol_lineage_order = c("VLMC", "OPC", "COP", "NFOL1", "NFOL2", "MFOL1",
        "MFOL2", "MOL1", "MOL2", "MOL3", "MOL4", "MOL5", "MOL6")
      tmp_df = tmp_df[match(ol_lineage_order, tmp_df$types) , ]
      tmp_df$types = factor(tmp_df$types, levels = unique(tmp_df$types))
      limits = aes(ymax = mean + se, ymin = mean - se)
      tmp_plot = ggplot(data = tmp_df, aes(x = types, y = mean)) + geom_point() +
        geom_errorbar(limits, width = 0.2) + theme_bw() +
        ylab("TMM-Normalized RNA CPM") + xlab("") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
          text = element_text(size = 22))
      print(tmp_plot)

    }

    get_marques_plot(input$gene)

  }

  if(input$data_type == "Tasic_RNA"){

    get_tasic_plot <- function(gene){

      if(!gene %in% rownames(tasic)){
        gene = cap_first(gene)
        if(!gene %in% rownames(tasic)){
          stop("Even when uncapitalized, that gene symbol is not in the Tasic et al. data set.")
        }
      }
      gene_data = tasic[gene, ]
      if(input$log){
        types = colnames(tasic)[grepl("log", colnames(tasic))]
        tmp = gene_data[ , types]
        mean = tmp[ , grepl("mean", colnames(tmp))]
        se = tmp[ , grepl("se", colnames(tmp))]
      } else {
        types = colnames(tasic)[!grepl("log", colnames(tasic))]
        tmp = gene_data[ , types]
        mean = tmp[ , grepl("mean", colnames(tmp))]
        se = tmp[ , grepl("se", colnames(tmp))]
      }
      types_clean = types[grepl("mean", types)]
      types_clean = sapply(strsplit(types_clean, "_", fixed = TRUE), "[", 1)
      tmp_df = data.frame(mean = as.numeric(mean), se = as.numeric(se), types = types_clean)
      tmp_df$types = factor(tmp_df$types, levels = unique(tmp_df$types))
      limits = aes(ymax = mean + se, ymin = mean - se)
      tmp_plot = ggplot(data = tmp_df, aes(x = types, y = mean)) + geom_point() +
        geom_errorbar(limits, width = 0.2) + theme_bw() +
        ylab("TMM-Normalized RNA CPM") + xlab("") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
          text = element_text(size = 22))
      print(tmp_plot)

    }

    get_tasic_plot(input$gene)

  }

  if(input$data_type == "Zeisel_RNA"){

    get_zeisel_plot <- function(gene){

      if(!gene %in% rownames(zeisel)){
        gene = cap_first(gene)
        if(!gene %in% rownames(zeisel)){
          stop("Even when uncapitalized, that gene symbol is not in the Zeisel et al. data set.")
        }
      }
      gene_data = zeisel[gene, ]
      if(input$log){
        types = colnames(zeisel)[grepl("log", colnames(zeisel))]
        tmp = gene_data[ , types]
        mean = tmp[ , grepl("mean", colnames(tmp))]
        se = tmp[ , grepl("se", colnames(tmp))]
      } else {
        types = colnames(zeisel)[!grepl("log", colnames(zeisel))]
        tmp = gene_data[ , types]
        mean = tmp[ , grepl("mean", colnames(tmp))]
        se = tmp[ , grepl("se", colnames(tmp))]
      }
      types_clean = types[grepl("mean", types)]
      types_clean = sapply(strsplit(types_clean, "_", fixed = TRUE), "[", 1)
      tmp_df = data.frame(mean = as.numeric(mean), se = as.numeric(se), types = types_clean)
      tmp_df$types = factor(tmp_df$types, levels = unique(tmp_df$types))
      limits = aes(ymax = mean + se, ymin = mean - se)
      tmp_plot = ggplot(data = tmp_df, aes(x = types, y = mean)) + geom_point() +
        geom_errorbar(limits, width = 0.2) + theme_bw() +
        ylab("TMM-Normalized RNA CPM") + xlab("") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1),
          text = element_text(size = 22))
      print(tmp_plot)

    }

    get_zeisel_plot(input$gene)

  }

  })

})
