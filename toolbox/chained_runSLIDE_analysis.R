#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(SLIDEHelper)
library(tidyverse)
library(doParallel)

yaml_path = args[1]


runSLIDEHelper_fromYamlpath = function(yaml_path) {

  yaml_input = yaml::yaml.load_file(yaml_path)

  er_path = yaml_input$out_path

  er_results_path = list.files(er_path, pattern = "final_delta_", full.names = T)

  if (length(er_results_path) == 0) {
    return()
  }
  er_results_temp = readRDS(er_results_path)

  x_path = yaml_input$x_path
  y_path = yaml_input$y_path

  x_mat = as.matrix(read.csv(x_path, row.names = 1, check.names = F))
  y_mat = as.matrix(read.csv(y_path, row.names = 1, check.names = F))

  z_mat = EssReg::predZ(x_mat, er_results_temp)
  colnames(z_mat) = paste0("Z", c(1:ncol(z_mat)))
  rownames(z_mat) = rownames(y_mat)
  z_path =  paste0(er_path, "/z_mat.csv")
  write.csv(z_mat, z_path)

  SLIDE_res <- runSLIDE(y_path = y_path,
                        z_path = z_path,
                        er_path = er_results_path,
                        do_interacts = TRUE,
                        spec = 0.1,
                        niter = 1000)

  num_top_feats <- 10
  condition <- "auc"
  SLIDE_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_results_path, out_path = er_path, SLIDE_res, num_top_feats = 10, condition)

  saveRDS(SLIDE_res, paste0(er_path, "/SLIDE_res.RDS"))

  plotSigGenes(SLIDE_res, out_path = er_path, plot_interaction = TRUE)

  #the SLIDE_res has to be the output from GetTopFeatures
  CalcControlPerformance(z_matrix = z_mat, y_path, SLIDE_res, niter = 1000, condition, out_path = er_path)

}



create_corr_network = function(x, y, out_dir, plot_name, sig_genes) {
  # change wd cause weird
  setwd(out_dir)

  x_gene = as.matrix(x[, sig_genes])

  col_auc = round(apply(x_gene, 2, function(xs) glmnet:::auc(as.matrix(y), as.matrix(xs))), 2)

  temp_cols = ifelse(col_auc > 0.55, "#40006D", ifelse(col_auc < 0.45, "#59A14F", "lightgray"))

  x_temp = cor(x_gene)

  pl = qgraph(x_temp, filename= plot_name,
              layout = "spring", threshold=0.4, repulsion=0.1,
              labels = stringr::str_trunc(colnames(x_temp), side = "right", width = 25), color = temp_cols,
              title = paste0("LF_", plot_name),
              label.scale.equal=FALSE,label.prop=0.95,shape="ellipse",
              posCol="salmon", negCol="skyblue",filetype='png',height=5,width=7)
}



create_corr_network_from_path = function(yaml_path) {
  yaml_input = yaml::yaml.load_file(yaml_path)

  er_path = yaml_input$out_path

  sig_genes_path = list.files(er_path, pattern = "plotSigGenes_data.RDS", full.names = T)

  if (length(sig_genes_path) == 0) {
    return()
  }
  sig_genes = readRDS(sig_genes_path)

  x_path = yaml_input$x_path
  y_path = yaml_input$y_path

  x_mat = as.matrix(read.csv(x_path, row.names = 1, check.names = F))
  y_mat = as.matrix(read.csv(y_path, row.names = 1, check.names = F))


  corr_dir = paste0(er_path, "/corr_networks")
  if (!dir.exists(corr_dir)) {
    dir.create(corr_dir)
  }

  for (lf in unique(sig_genes$lf_num)) {

    lf_data = sig_genes %>% filter(lf_num == lf)
    create_corr_network(x_mat, y_mat, corr_dir, lf, lf_data$names)
  }
}




runSLIDEHelper_fromYamlpath(yaml_path)
create_corr_network_from_path(yaml_path)
