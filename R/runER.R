#' @export
runER = function(yaml_path, coarseGrid = F, cleanData = T, removeCleanedData = T,
                lambdas = NULL, deltas = NULL, ...) {

  # check lambdas and deltas
  if (!is.null(lambdas) && any(lambdas) <= 0) {
    lambdas = c(1.0, 0.1)
  }
  
  if (!is.null(deltas) && any(deltas) <= 0) {
   deltas = c(0.1, 0.01) 
  }
  if (cleanData) {
    cleaned = JishnuLabTools::load_data(yaml = yaml_path, ...)
    yaml_path = cleaned$yaml
  }
  er_input = yaml::yaml.load_file(yaml_path)
  original_path = er_input$out_path

  if (coarseGrid) {
    er_input$k = JishnuLabTools::set_k_folds(er_input)
      
    if (is.null(deltas) && length(er_input$delta) == 1) {
      # default deltas and lambdas if we aren't passed list of deltas
      deltas = c(0.25, 0.1, 0.05, 0.01)
    } else if (is.null(deltas)) {
      # use the passed list of deltas
      deltas = er_input$delta
    }

    # check if we have a list of lambas (we should test more than 1, but dont
    # need to test many)
    if (is.null(lambdas) && length(er_input$lambda) > 1) {
      # set lambdas to this list
      lambdas = er_input$lambda
      # otherwise we just use the lambdas defined above
    }

    run_list = list()

    for (d in deltas) {
      for (l in lambdas) {
        er_input$out_path = paste0(original_path, 'lambda', l, "_delta", d, "/")

        run_list = append(run_list, er_input$out_path)

        er_input$delta = d
        er_input$lambda = l
        yaml::write_yaml(er_input, yaml_path)

        cat("running pipeline 3\n")

        EssReg::pipelineER3(yaml_path)

        # run slide on these results
        JishnuLabTools::run_slide(loaded_yaml = er_input)

        JishnuLabTools:::unregister_dopar()
        doParallel::registerDoParallel(cores)
      }
    }
  } else {
    yaml::write_yaml(er_input, yaml_path)

    EssReg::pipelineER3(yaml_path)
    JishnuLabTools::run_slide(loaded_yaml = er_input)
  }

  # add main direcotry and runs to run list
  run_list = list(run_dir = original_path, runs = run_list)
  # save run list
  saveRDS(run_list, paste0(original_path, "run_list.RDS"))

  if (cleanData == T && removeCleanedData == T) {
    if (file.exists(yaml_path$x_path) == T && file.exists(yaml_path$y_path) == T) {
    cat(" Deleting x & y saved in output directory\n")
    file.remove(yaml_path$x_path)
    file.remove(yaml_path$y_path)
    }
  }
}
