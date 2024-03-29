#' @export
run_slide = function(yaml_path = NULL, loaded_yaml = NULL, spec = 0.1, do_interacts = F) {

#   # check if we've registered parallel cores
#   running_parallel = exists(x = 'cores')
#   cores = ifelse(running_parallel, cores, NULL)


  if ( !is.null(loaded_yaml)) {
    yaml_input = loaded_yaml
  } else if ( !is.null(yaml_path)) {
    lower_yaml_path = stringr::str_to_lower(yaml_path)
    if (stringr::str_ends(lower_yaml_path, ".rds")) {
      # load as RDS
      yaml_input = readRDS(yaml_path)
    } else {
      yaml_input = yaml::yaml.load_file(yaml_path)
    }
  }

  results_folder = yaml_input$out_path

  results_folder_files = list.files(results_folder, full.names = T)

  # er results
  if(all(which(stringr::str_detect(results_folder_files, pattern = "final_"))) == F) {
    return()
  }
  er_results_path = results_folder_files[which(stringr::str_detect(results_folder_files, pattern = "final_delta_"))]

  er_results = readRDS(er_results_path)

  ######### run SLIDE

  x_mat = as.matrix(read.csv(yaml_input$x_path, row.names = 1))
  y_mat = as.matrix(read.csv(yaml_input$y_path, row.names = 1))


  # threshold A matrix
  new_A = EssReg::threshA(er_results$A, 0.08, FALSE)

  x_mat = scale(x_mat, T, T)

  er_results$A = new_A

  z_mat = EssReg::predZ(x_mat, er_results)
  colnames(z_mat) = paste0("Z", c(1:ncol(z_mat)))

  f_size_val = er_results$K
  if (f_size_val > 100) {
    f_size_val = 100
  } else if (f_size_val == 1) {
    cat("only one latent factor, skipping\n")
    return()
  }

  slide_res = SLIDE::SLIDE(z_mat, y_mat, method = 4, do_interacts = do_interacts, betas = NULL, top_prop = NULL, marginals = NULL,
                    spec = spec, fdr = 0.1, niter = 5000, elbow = FALSE, f_size = f_size_val, parallel = T, ncore = cores)

  saveRDS(slide_res, paste0(yaml_input$out_path, 'slide_res.RDS'))

  write.csv(z_mat, paste0(yaml_input$out_path, 'z_mat.csv'), row.names = T, col.names = T)

  marg_vars = slide_res$marginal_vars

  if (length(marg_vars) == 0) {
    cat("no significant latent factors, skipping\n")
    return()
  }

  ks = sapply(marg_vars, function(x) as.double(gsub("z", "", x, perl = T)))

  sig_genes = list()

  if (length(ks) > 0) {
    for (i in 1:length(ks)){
      feature_res <-SLIDE::SigGenes_A(er_results, ks[[i]], as.matrix(x_mat), y_mat, thre=0.1, negcol="blue", posCol="red",orderbycol=T)
      feature_res <- feature_res[order(feature_res[, 1], decreasing = TRUE), ]
      sig_genes[[i]] <- feature_res
    }
  }

  saveRDS(sig_genes, paste0(yaml_input$out_path, 'sig_genes.RDS'))


  # get significant genes plot
  sig_genes_res = JishnuLabTools:::plotSigGenes(sig_genes = sig_genes,
                                               er_input = er_results, slide_res = slide_res,
                                               xdf = x_mat, ydf = y_mat,
                                               output_plot_path = yaml_input$out_path)


  # ######### plot performance
  JishnuLabTools:::plot_performance(results_folder)
}
