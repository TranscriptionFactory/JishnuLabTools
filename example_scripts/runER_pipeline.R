#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(devtools)
library(doParallel)
library(foreach)
library(tidyverse)
library(optparse)

# if need to install
# main directories are Hanxi-002/EssReg and jishnu-lab/SLIDE should be installed here - 
devtools::install_github(repo = "Hanxi-00/EssReg", force = F, dependencies = T)
devtools::install_github(repo = "jishnu-lab/SLIDE", force = F, dependencies = T)

#devtools::install_github(repo = "TranscriptionFactory/JishnuLabTools", force = F, dependencies = T)

library(EssReg)
library(SLIDE)
#library(JishnuLabTools)

cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# if(!is.na(cores) & cores > 1) cores <- cores
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

# process arguments from command line

runER_args = list(
  optparse::make_option(c("-x", "--x_path"), action = "store", type = "character",
                        default = "", help = "Path to X data"),
  optparse::make_option(c("-y", "--y_path"), action = "store",type = "character",
                        default = "", help = "Path to Y data"),
  optparse::make_option(c("-d", "--data_path"), action = "store",type = "character",
                        default = "", help = "Path to data"),
  optparse::make_option(c("-p", "--yaml_path"), action = "store",type = "character",
                        default = "", help = "Path to yaml file"),
  optparse::make_option(c("-c", "--coarse_grid"), action = "store_true",
                        default = F, help = "Run coarse grid search of lamda/deltas"),
  optparse::make_option(c("-e", "--pipeline"), action = "store", type = "integer",
                        default = 3, help = "Pipeline # (1, 2, or 3)")
  )
command_args = optparse::parse_args(optparse::OptionParser(option_list = runER_args), args = args)

yaml_path = command_args$yaml_path

pipeline = command_args$pipeline

# run pipeline #

if (pipeline == 1) {
  # pipeline 1
  EssReg::pipelineER1(yaml_path)
} else if (pipeline == 2) {
  # pipeline 2
  EssReg::pipelineER2(yaml_path)
} else {
  # default is just to run 3
  # pipeline 3
  EssReg::pipelineER3(yaml_path)
}
