#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(devtools)
library(JishnuLabTools)

cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# if(!is.na(cores) & cores > 1) cores <- cores
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

getModelName = function(pth) {
  fname = tail(str_split(pth, pattern = "/")[[1]], 1)

  mdl = str_split(fname, pattern = ".rds")[[1]][1]
  return(mdl)
}

percentage_zeros = function(df, axs = 2, prop = 0.5) {

  prop_zeros_col = apply(df, axs, function(x) length(which(x == 0)) / length(x))

  return(which(prop_zeros_col < prop))
}

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
                        default = 3, help = "Pipeline # (1, 2, or 3)",
  optparse::make_option(c("-c", "--corr_filter"), action = "store", type = "double",
                        default = 0.0, help = "Filter correlated variables. Correlation threshold"),
  optparse::make_option(c("-o", "--out_path"), action = "store", type = "character",
                        default = ".", help = "Output path"),
  optparse::make_option(c("-m", "--filter_sparsity"), action = "store", type = "double",
                        default = 0.0, help = "Sparsity filter (0 = filter none, 0.5 = filter median, 1.0 = filter all"),
  optparse::make_option(c("-v", "--filter_var"), action = "store", type = "double",
                        default = 0.0, help = "Sparsity variance quantile (0 = filter none, 0.5 = filter median, 1.0 = filter all"))
)
command_args = optparse::parse_args(optparse::OptionParser(option_list = runER_args), args = args)



rds_path = command_args$data_path

# load yaml
yml = yaml::yaml.load_file(command_args$yaml_path)

out_path = command_args$out_path

cat("model ", getModelName(rds_path), "\n")
out_path = paste0(out_path, getModelName(rds_path), "/")

yml$out_path = out_path

cat("output path is ", out_path, "\n")

df = readRDS(rds_path)
rownames(df) = df[, 1]
df = df[, -1]

df = df[sample(1:nrow(df)), ]

# remove GM42418 and Lars2
bad_cols = which(stringr::str_detect(stringr::str_to_upper(names(df)), pattern = "GM42418|LARS2"))

if (length(bad_cols) > 0) {
  df = df[, -bad_cols]
}


# transcriptomics vs metabolomics

mode = ifelse(nrow(df) > 100, "TRX", "MTX")

cat("Dim df is pre-filter ", dim(df), "\n")

# replace all NAs with zeros
df[is.na(df)] = 0

# then remove columns with all zeros
df = df[, which(colSums(df) != 0)]

df_temp = df[, -1]

# remove cols with zero standard deviation
df_temp = df_temp[, which(apply(df_temp, 2, sd) != 0)]

if ( command_args$filter_sparsity > 0.0 ) {
  # need to filter out cells that have >50% zeros in each group
  non_sparse_rows = percentage_zeros(df_temp, axs = 1, prop = command_args$filter_sparsity)

  df_temp = df_temp[non_sparse_rows, ]

  # also need to filter y
  df = df[non_sparse_rows, ]

  # get percentage zeros
  non_sparse_cols = percentage_zeros(df_temp, axs = 2, prop = command_args$filter_sparsity)

  # keep only these
  df_temp = df_temp[, non_sparse_cols]
}

if ( command_args$filter_var > 0.0 ) {
  # filter by variance
  df_temp_col_var = apply(df_temp, 2, var)

  df_high_var = which(df_temp_col_var > quantile(df_temp_col_var, command_args$filter_var))

  df_temp = df_temp[, df_high_var]
}

# first filter by correlation
if ( command_args$corr_filter > 0.0 ){

  df_temp_cov = abs(cor(as.matrix(df_temp)))

  df_temp_cov[lower.tri(df_temp_cov, diag = T)] = 0

  df_cov_high = apply(df_temp_cov, 2, function(x) length(which(x > command_args$corr_filter)))

  # get quantile of # correlated
  df_temp = df_temp[, -which(df_cov_high > quantile(df_cov_high, 0.90))]
}

df = cbind.data.frame(df[, 1], df_temp)
cat("Dim df is post-filter ", dim(df), "\n")


# don't quantile filter columns if there aren't many
data = JishnuLabTools::load_data(obj = df, yaml = yml, remove_mito_ribo = F,
                                 quantile_filter = 0.0, remove_zero_median_cols = F)


er_input =  yaml::yaml.load_file(data$yaml)

if (nrow(df) < 20) {
  er_input$k = nrow(df)
} else if (nrow(df) < 50) {
  er_input$k = 5
} else {
  er_input$k = 10
}

yaml_path = data$yaml

if (nrow(df) > 100) {
  lbds = c(1.5, 1.0, 0.5)
  dlts = c(0.3, 0.2, 0.15, 0.1)
} else {
  lbds = c(0.75, 1.0)
  dlts = c(0.15, 0.2, 0.25)
}


orig_path = er_input$out_path

for (ind1 in 1:length(lbds)) {
  for (ind2 in 1:length(dlts)) {
    er_input$out_path = orig_path

    er_input$lambda = lbds[ind1]
    er_input$delta = dlts[ind2]

    er_input$out_path = stringr::str_c(er_input$out_path, "r_lambd", lbds[ind1], "_delt", dlts[ind2], "/", sep = "")
    yaml::write_yaml(er_input, yaml_path)

    cat("running pipeline 3 with reps \n")

    EssReg::pipelineER3(yaml_path)

    # run slide on these results
    JishnuLabTools::run_slide(loaded_yaml = er_input)

    yaml::write_yaml(er_input, yaml_path)

    unregister_dopar()
    registerDoParallel(cores)
  }
}




