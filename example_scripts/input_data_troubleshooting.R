#!/usr/bin/env Rscript

library(devtools)
# if need to install
# main directories are Hanxi-002/EssReg and jishnu-lab/SLIDE should be installed here -
devtools::install_github(repo = "Hanxi-002/EssReg", force = F, dependencies = T)
devtools::install_github(repo = "jishnu-lab/SLIDE", force = F, dependencies = T)

library(doParallel)
library(foreach)
library(tidyverse)
library(EssReg, attach.required = T)
library(SLIDE, attach.required = T)
library(Matrix)

###############################################################################
# Load Data
###############################################################################
# load your yaml file
yaml_input = yaml::yaml.load_file('path/to/you/yaml.yaml')

# load x and y
x = read.csv(yaml_input$x_path, row.names = 1)
y = read.csv(yaml_input$y_path, row.names = 1)


###############################################################################
# Filtering attempt 1: Try removing median = 0 features
###############################################################################
no_median_data = clean_data(xdata = x, ydata = y,
                            remove_zero_median_cols = T)

# save new x and y to folder
write.csv(no_median_data$x, 'path/to/where/you/want/to/save.csv')
write.csv(no_median_data$y, 'path/to/where/you/want/to/save.csv')


###############################################################################
# Filtering attempt 2: Try removing median = 0 features and low variance features
###############################################################################
# if this doesn't work, also remove the lowest 25% features by variance
high_var_data = clean_data(xdata = x, ydata = y,
                           remove_zero_median_cols = T,
                           col_var_quantile_filter = 0.25)

write.csv(high_var_data$x, 'path/to/where/you/want/to/save.csv')
write.csv(high_var_data$y, 'path/to/where/you/want/to/save.csv')



###############################################################################
# Update yaml file
###############################################################################
# update yaml and save
yaml_input$x_path = 'path/to/new/x/data'
yaml_input$y_path = 'path/to/new/y/data'

# save updated yaml

yaml::write_yaml(yaml_input, 'new/yaml/path.yaml')


###############################################################################
# Run ER
###############################################################################
# initialize cores for running ER
cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# if(!is.na(cores) & cores > 1) cores <- cores
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')


# run ER pipeline 3
EssReg::pipelineER3('new/yaml/path.yaml')








