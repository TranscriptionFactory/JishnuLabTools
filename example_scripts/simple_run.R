#!/usr/bin/env Rscript
args = R.utils::commandArgs(asValues = TRUE, excludeReserved = TRUE, args = TRUE)

# Arguments are:
#   -x path to x
#   -y path to y
#   -o output path
#   -e evaluation type (use 'r' for regression and 'c' for classication)
#   -c classification
#   -d path to data that is [Y, X] (e.g. Y in first column)
#   -a run all 3 pipelines
#   -p yaml parameters
#
#
#
#
arglist = list(x_path = list(options = c("x", "x_path"), default = ""),
               y_path = list(options = c("y", "y_path"), default = ""),
               out_path = list(options = c("o", "out_path"), default = "/"),
               run_all = list(options = c("a", "run_all"), default = F),
               yaml_file = list(options = c("p", "parameters", "yaml_parameters",
                                            "yaml_file", "yaml"), default = F),
               eval_type = list(options = c("r", "c", "regression",
                                            "classification", "auc", "corr"),
                                default = "regression"))


# always check for updates
devtools::install_github("TranscriptionFactory/JishnuLabTools",
                         force = F,
                         dependencies = TRUE)

devtools::install_github("Hanxi-002/EssReg",
                         force = F,
                         dependencies = TRUE)

library(JishnuLabTools)
library(EssReg)


# load arguments
loaded_args = JishnuLabTools::arg_loader(arglist)

if (loaded_args$yaml_file == F) {
  # load yaml
  if (loaded_args$eval_type %in% c("r", "corr", "regression")) {
    yaml = JishnuLabTools::regression
  } else if (loaded_args$eval_type %in% c("c", "auc", "classification")) {
    yaml = JishnuLabTools::classification
  } else {
    cat("Error with eval_type argument. Pass yaml file or check arguments\n")
  }
} else {
  # load the passed yaml
  yaml = yaml::yaml.load_file(loaded_args$yaml_file)
}




cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# if(!is.na(cores) & cores > 1) cores <- cores
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')



unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

