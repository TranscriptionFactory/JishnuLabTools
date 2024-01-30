#!/usr/bin/env Rscript

# '''
#   instructions:
#     1. save a dummy yaml file in path w/ correct parameters for classification (eval: auc, levels: [0,1]) vs regression (eval: corr, levels: NULL)
#       - search for 'DUMMY YAML FILE'
#     2. search for 'CLUSTERING TUNABLE PARAMETERS' and change as necessary
#     3. search for 'SLIDE TUNABLE PARAMETERS' and change as necessary
#     4. enjoy :)
# '''


# this script runs ER without cross validation and instead runs cross validation during SLIDE; much faster & less error prone
# usage w/ slurm:
# ''' Rscript path/to/this/script.R path/to/RDS/with/Y_labels/in_first_column_and_features_after/ path/to/yaml/file.yaml '''

# install most recent packages like this:
# '''
# devtools::install_github("jishnu-lab/SLIDE")
# devtools::install_github("Hanxi-002/SLIDEHelper")
# '''

#####################################
# begin: CLUSTERING TUNABLE PARAMETERS
# These parameters are important - make sure to check before running; they are fairly conservative so you should be fine
num_k_folds = 5 # change as necessary
threshold_fdr = 0.05 # false discovery rate threshold
alpha_level = 0.1 # alpha level for significance
yaml_path = '.' # REPLACE ME WITH PATH TO DUMMY YAML FILE
scale_data_mean_zero_std_dev_1 = F # make true if you want to scale data to have mean zero and std dev = 1
# end: CLUSTERING TUNABLE PARAMETERS
#####################################

# this just names the model with the date/time the script was run; feel free to change; replace spaces just in case they
# want to cause problems w/ file naming later
model_name = stringr::str_replace_all(Sys.time(), pattern = " ", replacement = "_")

args = commandArgs(trailingOnly=TRUE) # parses command line arguments

# or,
rds_path = args[1] # replace rds_path w/ path to data (y = first column, x = features in following columns)
out_path = paste0(args[2], "/") # and out_path (where you want results saved)

# a fairly comprehensive array of lambda/delta parameters (will do # lambdas * # delta runs)
lbds = c(1.0, 0.5, 0.1)
dlts = c(0.3, 0.2, 0.1, 0.2, 0.05, 0.01)

library(matrixcalc)
library(ROCR)
library(e1071)
library(doParallel)
library(foreach)
library(scales)
library(doRNG)
library(matlib)
library(BAS)
library(gmp)
library(tidyverse)
library(devtools)
library(EssReg)
library(ggpubr)
library(SLIDE)
library(SLIDEHelper)

cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

cat("Path to RDS file is ", rds_path, "\n")

df = readRDS(rds_path)

out_path = paste0(out_path, "/", model_name, "/")

# load yaml
yml = yaml::yaml.load_file(yaml_path)

yml$out_path = out_path


if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = T)
}


x_path = paste0(out_path, "x.csv")
y_path = paste0(out_path, "y.csv")

x_data = df[, -1]

# uncomment below if you want to scale each column to have mean = 0 and std dev = 1
if (scale_data_mean_zero_std_dev_1) {
  x_data = apply(x_data, 2, function(x) scale(x, T, F))
}

write.csv(df[, 1], y_path)
write.csv(x_data, x_path)


cat('Saving X and Y to path: ', out_path, "\n")

yml$x_path = x_path
yml$y_path = y_path

yml$k = num_k_folds
yml$thresh_fdr = threshold_fdr
yml$alpha_level = alpha_level
yml$lasso = F

er_input = yml

orig_path = paste0(er_input$out_path, "pipeline3/")
yaml_path = paste0(out_path, "params.yaml")
# save yaml
yaml::write_yaml(yml, yaml_path)

for (ind1 in 1:length(lbds)) {
  for (ind2 in 1:length(dlts)) {
    er_input$out_path = orig_path

    er_input$lambda = lbds[ind1]
    er_input$delta = dlts[ind2]

    er_input$out_path = stringr::str_c(er_input$out_path, "r_lambd", lbds[ind1], "_delt", dlts[ind2], "/", sep = "")
    yaml::write_yaml(er_input, yaml_path)
    cat("running pipeline 3 with reps \n")

    EssReg::parseRun(yaml_path)

    # RUN SLIDE
    ###################################
    yaml_input = yaml::yaml.load_file(yaml_path)

    er_path = yaml_input$out_path

    er_results_path = list.files(er_path, pattern = "final_", full.names = T)

    if (length(er_results_path) == 0) {
      cat("\nno er results found, exiting\n")
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

#####################################
# begin: SLIDE TUNABLE PARAMETERS

    # add all the stuff to our orig yaml that we'll need for cv
    er_input$z_path = z_path
    er_input$niter = 20
    er_input$parallel = TRUE
    er_input$method = 4
    er_input$ncore = cores
    er_input$fdr = threshold_fdr
    er_input$spec = 0.1
    er_input$f_size = 20

# end: SLIDE TUNABLE PARAMETERS
#####################################

    yaml::write_yaml(er_input, yaml_path)

    cat("beginning SLIDE \n")

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
    ###################################
    cat("beginning SLIDE CV \n")

    # SLIDE CV
    ##################################


    for(j in 1:er_input$nreps){

      paperBench(yaml_path, replicate=j)

    }

    er_input <- yaml::yaml.load_file(yaml_path)

    pathLists <- list.files(er_input$out_path,recursive = T,pattern = "results")
    perfList <- lapply(paste0(er_input$out_path,pathLists), readRDS)
    perRes <- do.call(rbind,lapply(perfList,function(x){x$final_corr}))


    if (er_input$eval_type == "corr") {

      lambda_boxplot = ggpubr::ggboxplot(data = perRes, x = "method", y = "corr", palette = "aaas",
                                         fill = "method" ) +
        ggpubr::stat_compare_means(label = "p.signif")

    } else {

      lambda_boxplot = ggpubr::ggboxplot(data = perRes, x = "method", y = "auc", palette = "aaas",
                                         fill = "method" ) +
        ggpubr::stat_compare_means(label = "p.signif")
    }

    ggplot2::ggsave(plot = lambda_boxplot, filename = paste0(er_input$out_path, "cv_boxplot.png"), height = 6, width = 6)

    saveRDS(perRes,file=paste0(er_input$out_path,"boxplot_data.rds"))

  }
}
