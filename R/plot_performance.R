#' @export
plot_performance = function(path_to_results, save_plot = T) {

  results_folder_files = list.files(path_to_results, full.names = T)

  er_model_eval_path = results_folder_files[which(str_detect(results_folder_files, pattern = "pipeline_step5"))]

  if ( length(er_model_eval_path) == 0 ){
    cat("No ER performance data to plot \n")
    return()
  }

  er_model_eval = readRDS(er_model_eval_path)

  er_model_eval$method = factor(er_model_eval$method, levels = c("plainER", "plainER_y",
                                                                 "lasso", "lasso_y"))
  evaltype = "auc"
  if (names(er_model_eval)[2] == "corr") {
    evaltype = "corr"
  }
  boxplot = ggpubr::ggboxplot(er_model_eval, x = "method", y = evaltype, palette = "npg",
                              fill = "method" ) +
    ggpubr::stat_compare_means(comparisons = list(c("plainER", "plainER_y"),
                                                  c("lasso", "lasso_y")))

  if ( save_plot ) {
    ggplot2::ggsave(plot = boxplot, filename = paste0(path_to_results, 'auc_boxplot.png'))
  }
  return( boxplot )
}
