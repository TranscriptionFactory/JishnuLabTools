#!/usr/bin/env Rscript
library(SLIDEHelper)
library(tidyverse)
library(doParallel)

# devtools::install_github("Hanxi-002/SLIDEHelper", force = F, dependencies = F)

cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')


runSLIDEHelper = function(yaml_path) {

  yaml_input = yaml::yaml.load_file(yaml_path)

  er_path = yaml_input$out_path

  er_results_path = list.files(er_path, pattern = "final_delta_", full.names = T)
  er_results_temp = readRDS(er_results_path)

  x_path = yaml_input$x_path
  y_path = yaml_input$y_path

  x_mat = as.matrix(read.csv(x_path, row.names = 1, check.names = F))
  y_mat = as.matrix(read.csv(y_path, row.names = 1, check.names = F))

  z_mat = EssReg::predZ(x_mat, er_results_temp)
  colnames(z_mat) = paste0("Z", c(1:ncol(z_mat)))
  rownames(z_mat) = rownames(y_mat)
  write.csv(z_mat, paste0(er_path, "/z_mat.csv"))

  SLIDE_res <- runSLIDE(y_path = y_path,
                        z_path = NULL,
                        z_matrix = z_mat,
                        er_path = er_results_path,
                        do_interacts = TRUE,
                        spec = 0.2,
                        niter = 1000)

  num_top_feats <- 10
  condition <- "auc"
  SLIDE_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_results_path, out_path = er_path, SLIDE_res, num_top_feats = 10, condition)

  saveRDS(SLIDE_res, paste0(er_path, "/SLIDE_res.RDS"))

  # get red/blue text plots
  plotSigGenes(SLIDE_res, out_path = er_path, plot_interaction = TRUE)

  #the SLIDE_res has to be the output from GetTopFeatures
  CalcControlPerformance(z_matrix = z_mat, y_path, SLIDE_res, niter = 1000, condition, out_path = er_path)

}
