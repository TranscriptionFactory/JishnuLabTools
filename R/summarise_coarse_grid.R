#' @export
summarise_coarse_grid = function(path) {

  default_name_order = T
  # get runs - try to find files ordered lambda_delta
  all_runs = grep(x = list.dirs(path),
                  pattern = "lambd(a)?([0-9]+)(\\.)?([0-9]*)_delt(a)?([0-9]+)(\\.)?([0-9]*)$",
                  value = T, ignore.case = T)

  if (length(all_runs) == 0) {
    default_name_order = F
    # maybe its delta_lambda (other order)
    all_runs = grep(x = list.dirs(path),
                    pattern = "delt(a)?([0-9]+)(\\.)?([0-9]*)_lambd(a)?([0-9]+)(\\.)?([0-9]*)$",
                    value = T, ignore.case = T)
  }


  runs_with_filename = grep(x = list.files(all_runs, full.names = TRUE,recursive = TRUE),
                            pattern = "plotSigGenes_data.RDS$",
                            value = T, ignore.case = T)

  valid_runs = stringr::str_split_i(runs_with_filename, i = 1,
                                    pattern = "plotSigGenes_data.RDS$")

  for (f in valid_runs) {
    res = data.frame()
    for (r in list.dirs(f, full.names = T, recursive = F)) {

      run_files = list.files(r, full.names = T)

      output_files = c("plotSigGenes_data.RDS", "pipeline_step5")

      plot_data = run_files[stringr::str_which(run_files, pattern = output_files[1])]

      model_eval_file = run_files[stringr::str_which(run_files, pattern = output_files[2])]
      # want red/blue balance to be higher than 90/10

      color_bal_blue = 0
      color_bal_red = 0
      num_lfs = 0
      total_genes = 0
      mean_auc = 0
      delta = 0
      lambda = 0

      if (length(plot_data) != 0) {
        # check sig genes
        plot_data = readRDS(plot_data)

        num_lfs = length(unique(plot_data$fac))

        color_bal_blue = length(which(plot_data$color == "Blue")) / length(plot_data$color)

        color_bal_red = length(which(plot_data$color == "Red")) / length(plot_data$color)

        total_genes = length(plot_data$color)
      }
      if (length(model_eval_file) != 0) {
        # get mean AUC
        mean_auc = 0
        model_eval_data = readRDS(model_eval_file)
        mean_auc = mean((model_eval_data %>% filter(method == "plainER"))$auc)


        split_run_name = stringr::str_split(r, pattern = "delt(a)?|lambd(a)?|_")[[1]]
        split_run_name = tail(split_run_name, 3)

        if (default_name_order) {
          lambda = split_run_name[1]
          delta = split_run_name[3]
        } else {
          lambda = split_run_name[3]
          delta = split_run_name[1]
        }


        res = rbind(res, data.frame(run = r,
                                    mean_auc = mean_auc,
                                    delta = delta,
                                    lambda = lambda,
                                    color_bal_blue = color_bal_blue,
                                    color_bal_red = color_bal_red,
                                    total_lf_genes = total_genes,
                                    num_lfs = num_lfs))
      }
    }
    saveRDS(res, paste0(f, "/summary.RDS"))

    if (length(res) > 0) {
      pl = ggplot(res, aes(lambda, delta, size = mean_auc, color = color_bal_red)) +
        geom_count() +
        scale_color_gradient(low = "blue", high = "red")

      pl_size = ggplot(res, aes(lambda, delta, size = (total_lf_genes)/num_lfs, color = num_lfs)) +
        geom_count() +
        scale_color_distiller(palette = "Spectral")

      combined_plot = ggpubr::ggarrange(plotlist = list(pl, pl_size), ncol = 2, nrow = 1)
      ggsave(plot = combined_plot, filename = paste0(f, "/combined_delta_lambda_grid.png"),
             device = "png", width = 12, height = 5, units = "in")
    }
  }
}
