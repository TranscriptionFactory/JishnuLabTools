
#' @export
smartER = function(yaml_path, all = F, run_repeats = F) {

  ################################################################################
  # load
  ################################################################################

  # load yaml file so we can output in slurm file and use later ifneeded
  er_input <- yaml::yaml.load_file(yaml_path)

  cat("YAML FILE")
  cat(yaml_path)
  cat("\n")

  # write yaml parameters to file
  for (listElement in 1:length(er_input)) {
    cat(names(er_input)[listElement])
    cat( " = ")
    cat(er_input[[listElement]])
    cat("\n")
  }


  getBestDeltaIndex = function(results) {
    # loop through results and find delta with highest value for plainER_true - plainER_permuted
    bestDelta = 0
    bestDeltaIndex = 0
    for (res in 1:length(results)) {
      if (!is.null(results[[res]]$result$auc)) {
        aucs = results[[res]]$result$auc
        grp_res = data.frame(results[[res]]$result) %>% group_by(method) %>% summarise(auc = mean(auc))

      } else if (!is.null(results[[res]]$result$corr)) {
        corrs = results[[res]]$result$corr
        grp_res = data.frame(results[[res]]$result) %>% group_by(method) %>% summarise(auc = mean(corr))

      }

      # just taking highest auc/corr
      deltaDiff = grp_res[3, 2]
      if (deltaDiff > bestDelta) {
        bestDelta = deltaDiff
        bestDeltaIndex = res
      }
    }
    bestDeltaIndex = ifelse(bestDeltaIndex > 0, bestDeltaIndex, 1)
    return(bestDeltaIndex)
  }

  ################################################################################
  # running ER
  ################################################################################

  new_yaml_path = stringr::str_replace(yaml_path, fixed(".yaml"), "_temp.yaml")

  original_path = er_input$out_path

  # create dir if doesn't exist
  if ( !dir.exists(er_input$out_path )) {
    dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)
  }

  if (all) {
    ###############################################################################
    # pipeline 1
    ###############################################################################

    cat("running pipeline 1\n")

    # change outpath to put results in pipeline1 directory

    er_input$out_path = paste0(original_path, 'pipeline1/')
    yaml::write_yaml(er_input, new_yaml_path)

    EssReg::pipelineER1(new_yaml_path, steps = "all")

    # read RDS from pipeline1 to pick approx best delta
    pipeline1_results = readRDS(paste0(er_input$out_path, 'pipeline_step2.rds'))

    bestDelta = getBestDeltaIndex(pipeline1_results)

    # use this delta for the next run
    nextDelta = readRDS(paste0(er_input$out_path, 'pipeline_step1.rds'))[[bestDelta]]$opt_delta

    # set new delta
    er_input$delta = nextDelta

    ###############################################################################
    # pipeline 2
    ###############################################################################

    # change outpath to put results in pipeline2 directory
    er_input$out_path = paste0(original_path, 'pipeline2/')

    # set lambdas
    er_input$lambda = c(0.1, 0.5, 1.0, 1.5)

    yaml::write_yaml(er_input, new_yaml_path)

    cat("running pipeline 2\n")

    # run pipeline2
    EssReg::pipelineER2(new_yaml_path, steps = "all")

    # find best delta and lambda from results
    # get optimal lambda and delta values from this index in step3 results
    nextDelta = readRDS(paste0(er_input$out_path, 'pipeline_step3.rds'))$opt_delta

    pipeline2_results = readRDS(paste0(er_input$out_path, 'pipeline_step4.rds'))

    # nextLambda = getBestDeltaIndex(pipeline2_results)

    nextLambda = readRDS(paste0(er_input$out_path, 'pipeline_step3.rds'))$opt_lambda
    # get opt lambda from cv

    # use these in pipeline3

    # set new delta and lambda
    er_input$delta = nextDelta
    er_input$lambda = er_input$lambda[nextLambda]
  }

  ###############################################################################
  # pipeline 3
  ###############################################################################

  # change outpath to put results in pipeline2 directory
  er_input$out_path = original_path
  yaml::write_yaml(er_input, new_yaml_path)

  if (runRepeats) {
    JishnuLabTools::runER(new_yaml_path, run_repeats = T)
  } else {
    JishnuLabTools::runER(new_yaml_path, run_repeats = F)
  }
  cat("Finished successfully\n")
}
