
library(SLIDEHelper)
library(EssReg)
library(tidyverse)
library(doParallel)

cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

cross_predict = function(er_results, orig_x, orig_y, orig_z,
                         new_x, new_y, slide_res = NULL) {

  # use original A matrix
  A_mat = er_results$A
  beta = er_results$beta

  if ( !is.null(slide_res) ) {

    sig_lfs = stringr::str_to_upper(slide_res$marginal_vars)

    non_sig_lf_cols = which(!colnames(A_mat) %in% sig_lfs)

    beta[non_sig_lf_cols] = 0

  } else {
    # sig_lfs = as.numeric(stringr::str_replace(slide_res$marginal_vars,
    #                                           pattern = "z", replacement = "")))

    sig_lfs = colnames(A_mat)
  }

  # get the names from x that are in new x
  orig_in_new = colnames(orig_x)[which(colnames(orig_x) %in% colnames(new_x))]
  orig_not_in_new = colnames(orig_x)[which(!colnames(orig_x) %in% colnames(new_x))]

  all_names_intersection = base::intersect(colnames(orig_x), colnames(new_x))

  new_x_worg = new_x


  cat("\n Total number of genes needed to make latent factors are n = ", length(colnames(orig_x)), " genes")
  if ( length(orig_in_new) > 0 ) {
    new_x_worig = new_x[, orig_in_new]
    cat("\n Using n = ", length(orig_in_new), " overlapping genes to predict with \n")
  } else {
    cat("\n No overlapping to predict with \n")
    return()
  }

  if ( length(orig_not_in_new) > 0 ) {
    to_add_to_x = matrix(data = 0, nrow = nrow(new_x_worig), ncol = length(orig_not_in_new))
    colnames(to_add_to_x) = orig_not_in_new
    new_x_worg = cbind(new_x_worig, to_add_to_x)

    # also make zero in original x
    orig_x[, orig_not_in_new] = 0

    cat("\n Adding n = ", length(orig_not_in_new), " columns of zeros for missing genes to make new X proper size \n")
  } else {
    new_x_worg = new_x_worig
  }

  reorder_y_to_xrows = function(x, y) {

    if (all(rownames(x) %in% rownames(y))) {
      reordered_y = as.matrix(y)[rownames(x), ]
      # return(as.matrix(reordered_y))
      return(reordered_y)
    }
    return(y)
  }


  # reorder colnames
  # new_x_worg = as.matrix(new_x_worg[, match(colnames(new_x_worg), colnames(orig_x))])
  new_x_worg = as.matrix(new_x_worg)

  # make sure new_y is properly ordered

  new_y = reorder_y_to_xrows(new_x_worg, new_y)

  # change loadings to just have sig Zs  #
  # # new_er_res$A[!c(colnames(new_er_res$A)) %in% sig_lfs] = 0
  # er_results = new_er_res

  # get new predicted Z
  new_z = EssReg::predZ(new_x_worg, er_results)

  # get new prediction
  new_y_pred = new_z %*% beta

  # remake z matrix so we use the same allocation for both
  # orig_y_pred = as.matrix(orig_z) %*% as.matrix(beta)
  orig_z_remade = EssReg::predZ(as.matrix(orig_x), er_results)

  orig_y_pred = as.matrix(orig_z_remade) %*% as.matrix(beta)

  is_classification = ifelse(length(unique(orig_y[, 1])) == 2, T, F)

  if (is_classification) {
    pred_auc = glmnet:::auc(as.matrix(new_y), new_y_pred)

    new_pred = ROCR::prediction(new_y_pred, new_y)
    new_perf = ROCR::performance(new_pred, "tpr", "fpr")

    # use orig Z and orig beta to get orig Y
    orig_pred_auc = glmnet:::auc(as.matrix(orig_y), orig_y_pred)

    orig_pred = ROCR::prediction(orig_y_pred, labels = orig_y)
    orig_perf = ROCR::performance(orig_pred, "tpr", "fpr")

    return(list(orig_perf = orig_perf, new_perf = new_perf,
                orig_eval = orig_pred_auc,
                new_eval = pred_auc))
  } else {

    orig_df = cbind.data.frame(orig_y, orig_y_pred)

    names(orig_df) = c("true", "pred")

    new_df = cbind.data.frame(new_y, new_y_pred)
    names(new_df) = c("true", "pred")

    return(list(orig_perf = orig_df, new_perf = new_df,
                orig_eval = cor(orig_y, orig_y_pred, method = "spearman"),
                new_eval = cor(new_y, new_y_pred)))
  }
}

cor_ggpubr = function(orig_perf, new_perf, title = NULL, orig_new_name = NULL) {

  if (is.null(orig_new_name)) {
    orig_new_name = c("ER", "crosspred")
  }

  # combine dfs
  orig_perf$run = orig_new_name[1]
  new_perf$run = orig_new_name[2]


  df = rbind.data.frame(orig_perf, new_perf)

  title = ifelse(!is.null(title), title, "ER Cross Prediction")

  pl = ggpubr:: ggscatter(df, x = "true",
                          y = "pred",
                          palette = "npg",
                          xlab = "True Y",
                          ylab = "Predicted Y",
                          facet.by = "run",
                          nrow = 2,
                          scales = "free",
                          add = "reg.line",
                          title = title,
                          add.params = list(color = "blue", fill = "lightgray"),
                          conf.int = T) +
    stat_cor(method = "spearman")

  return(pl)
}


perf_ggplot = function(orig_perf, new_perf, title = NULL,
                       orig_new_name = NULL) {


  if (is.null(orig_new_name)) {
    orig_new_name = c("ER", "crosspred")
  }

  make_df_from_ROCR_performance = function(perf) {

    return(list2DF(list(x = unlist(perf@x.values), y = unlist(perf@y.values))))
  }

  # make into dataframe
  orig_df = make_df_from_ROCR_performance(orig_perf)
  new_df = make_df_from_ROCR_performance(new_perf)

  # add annotation column for source of prediction
  orig_df$run = orig_new_name[1]
  new_df$run = orig_new_name[2]

  df = rbind(orig_df, new_df)

  title = ifelse(!is.null(title), title, "ER Cross Prediction")

  pl = ggplot2::ggplot(data = df, aes(x = x, y = y, color = run)) +
    ggplot2::geom_step(linewidth = 2) + theme_bw(base_size = 12) +
    xlab(orig_perf@x.name) + ylab(orig_perf@y.name) + ggtitle(title)


  df$x = round(df$x, 1)
  pl = ggpubr::ggline(data = df, x = "x", y = "y", color = "run", plot = "l",
                      size = 3.25, ylab = "True Positive Rate",
                      xlab = "False Positive Rate",
                      palette = "aaas", title = title)



  return(pl)
}

recluster_cross_pred = function(orig_run_path, new_run_path,
                                out_path = NULL,
                                orig_new_name = NULL, delta = 0.1,
                                lambda = 1.0, spec = 0.1, ...) {

  args = list(...)

  # load x,y we want to use
  orig_x = read.csv(paste0(orig_run_path, "/x.csv"), row.names = 1)
  orig_y = read.csv(paste0(orig_run_path, "/y.csv"), row.names = 1)

  # load new X, Y
  new_x = read.csv(paste0(new_run_path, "/x.csv"), row.names = 1)
  new_y = read.csv(paste0(new_run_path, "/y.csv"), row.names = 1)


  # set feature names
  colnames(orig_x) = stringr::str_to_upper(colnames(orig_x))
  colnames(new_x) = stringr::str_to_upper(colnames(new_x))

  all_names_intersection = base::intersect(colnames(orig_x),  colnames(new_x))

  run_data = list()

  orig_x = orig_x[, all_names_intersection]
  new_x = new_x[, all_names_intersection]

  run_data$orig_x = orig_x
  run_data$orig_y = orig_y
  run_data$new_x = new_x
  run_data$new_y = new_y

  is_classification = ifelse(length(unique(orig_y[, 1])) == 2, T, F)

  # run ER on data we want to use
  er_results = do.call(EssReg::plainER(x = orig_x,
                                       y = orig_y,
                                       sigma = NULL,
                                       delta = delta,
                                       lambda = lambda), args)


  er_results_path = paste0(out_path, "/orig_er_res.RDS")
  saveRDS(er_results, er_results_path)

  run_data$er_results = er_results

  saveRDS(run_data, paste0(out_path, "/input_data.RDS"))

  orig_z = EssReg::predZ(orig_x, er_results)
  colnames(orig_z) = paste0("Z", c(1:ncol(orig_z)))
  rownames(orig_z) = rownames(y_mat)
  z_path =  paste0(out_path, "/orig_z.csv")
  write.csv(orig_z, z_path)

  # write y
  y_path = paste0(out_path, "/orig_y.csv")
  write.csv(orig_y, y_path)

  x_path = paste0(out_path, "/orig_x.csv")
  write.csv(orig_x, x_path)

  SLIDE_res <- do.call(SLIDEHelper::runSLIDE(y_path = y_path,
                                             z_path = z_path,
                                             er_path = er_results_path,
                                             do_interacts = TRUE,
                                             spec = spec,
                                             niter = niter), args)


  if (is_classification) {
    condition <- "auc"
  } else {
    condition = "corr"
  }

  SLIDE_res <- SLIDEHelper::GetTopFeatures(x_path, y_path, er_results_path, out_path = out_path,
                                           SLIDE_res, num_top_feats = 10, condition)

  saveRDS(SLIDE_res, paste0(out_path, "/SLIDE_res.RDS"))

  SLIDEHelper::plotSigGenes(SLIDE_res, out_path = out_path, plot_interaction = TRUE)

  #the SLIDE_res has to be the output from GetTopFeatures
  SLIDEHelper::CalcControlPerformance(z_matrix = z_mat, y_path, SLIDE_res, niter = 1000,
                                      condition, out_path = out_path)


  plot_title = ifelse(is_classification, "Generalizing LFs: AUC for ",
                      "Generalizing LFs: Corr for ")

  # load SLIDE results if we want to cross predict using only significant LFs
  if (load_slide) {
    slide_res = SLIDE_res

    plot_title = paste0(plot_title, "SLIDE\n")
  } else {
    slide_res = NULL
    plot_title = paste0(plot_title, "ER\n")

  }

  perf_data = cross_predict(orig_er, orig_x, orig_y, orig_z, new_x, new_y, slide_res)

  perf_data$orig_eval = round(perf_data$orig_eval, 2)
  perf_data$new_eval = round(perf_data$new_eval, 2)

  # make a title with specific names if we have them

  if (!is.null(orig_new_name) & length(orig_new_name) == 2) {
    # replace title with names for groups
    plot_title = paste0(plot_title, orig_new_name[1], ": ",
                        perf_data$orig_eval,
                        "\n", orig_new_name[2], ": ", perf_data$new_eval)
  } else {
    plot_title = paste0(plot_title, "ER: ",
                        perf_data$orig_pred_auc, "\nCross Prediction: ",
                        perf_data$new_pred_auc)
  }


  # call plotting functions
  if (is_classification) {
    pl = perf_ggplot(perf_data$orig_perf, perf_data$new_perf, title = plot_title,
                     orig_new_name = orig_new_name)
  } else {

    pl = cor_ggpubr(perf_data$orig_perf, perf_data$new_perf, title = plot_title,
                    orig_new_name = orig_new_name)
  }

  save_path = ifelse(!is.null(out_path), out_path, new_run_path)

  # prepend name based on whether we used slide results
  save_path = ifelse(load_slide, paste0(save_path, "/sigLFs_"),
                     paste0(save_path, "/allLFs_"))

  # save plot and data
  saveRDS(perf_data, paste0(save_path, "cross_pred_performance.RDS"))
  ggplot2::ggsave(filename = paste0(save_path, "crosspred_perf_plot.pdf"), plot = pl,
                  height = 6, width = 5.75)
}






