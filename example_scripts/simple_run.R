#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(devtools)
library(doParallel)
library(foreach)
library(tidyverse)

# if need to install
devtools::install_github(repo = "TranscriptionFactory/JishnuLabTools", force = F, dependencies = T)

library(JishnuLabTools)

cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# if(!is.na(cores) & cores > 1) cores <- cores
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')


# process arguments from command line
command_args = optparse::parse_args(optparse::OptionParser(option_list = JishnuLabTools::runER_args), args = args)

yaml_path = command_args$yaml_path
coarseGrid = command_args$coarse_grid

# call ER function
JishnuLabTools::runER(yaml_path, coarseGrid)
