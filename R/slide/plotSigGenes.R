#' @export
plotSigGenes = function(sig_genes = NULL, er_input = NULL,
                        slide_res = NULL, xdf = NULL, ydf = NULL) {


  sdf = JishnuLabTools::check_for_df_or_path(sig_genes)

  er_results = JishnuLabTools::check_for_df_or_path(er_input)

  slide_results = JishnuLabTools::check_for_df_or_path(slide_res)

  x = JishnuLabTools::check_for_df_or_path(xdf)

  y = JishnuLabTools::check_for_df_or_path(ydf)

  if(any(is.null(c(sdf, x, y)))) {
    cat("Error loading data")
    return()
  }

  ks = as.numeric(stringr::str_remove(slide_results$marginal_vars, "z"))

  max_num_to_plot = 1

  cdf = data.frame()

  # test top genes
  sdf = JishnuLabTools::getTopGenes(er_results, ks, x, y)

  for (s in 1:length(sdf)) {
    temp = sdf[[s]]

    num_to_plot = ifelse(nrow(temp) < 20, nrow(temp), 20)
    max_num_to_plot = ifelse(num_to_plot > max_num_to_plot, num_to_plot, max_num_to_plot)

    temp$heights = seq(1, num_to_plot)
    temp$fac = s

    cdf = rbind.data.frame(cdf, temp)
  }

  plt = cdf %>% ggplot(., aes(x = factor(fac), y = heights, label = names)) +
    geom_text(aes(color = factor(color))) +
    scale_color_manual(values = c("blue", "red"), guide = "none") + theme_void() +
    theme(axis.text.x = element_text(), axis.title.x = element_text(),
          axis.title.y = element_text(angle = 90)) +
    xlab("Latent Factor") + ylab("Genes Associated with Latent Factor Cluster") + ylim(0, 20) +
    ggtitle("Genes Associated with Latent Factor")
  return(list("plt" = plt, "plot_df" = cdf))
}
