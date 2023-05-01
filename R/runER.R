#' @export
runER = function(yaml_path, runRepeats = F) {

  er_input = yaml::yaml.load_file(yaml_path)
  original_path = er_input$out_path

  er_input$out_path = paste0(original_path, 'pipeline3/')

  if (runRepeats) {
    er_input$k = 5
    lambdas = c(1.0, 0.1)
    deltas = c(0.25, 0.1, 0.05, 0.01)

    for (l in lambdas) {
      for (d in deltas) {
        er_input$out_path = paste0(original_path, 'pipeline3/delta', d, "_lambda", l, "/")
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
}
