library(shiny)
library(ggplot2)

shinyServer(function(input, output) {

  output$cellPlot <- renderPlot(height = 600, units = "px", width = 600, {

  std_error <- function(x) sd(x)/sqrt(length(x))

  if(input$data_type == "Sharma_RNA"){

    get_sharma_rpkm_plot <- function(gene){
      if(!gene %in% rownames(sharma_rpkm)){
        stop("That gene symbol is not in the Sharma et al. RPKM data set.")
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
        stop("That gene symbol is not in the Sharma et al. LFQ data set.")
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

  })

})
