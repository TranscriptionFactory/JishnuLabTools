#' @export
plotSigGenes = function(sig_genes = NULL, er_input = NULL,
                        slide_res = NULL, xdf = NULL, ydf = NULL,
                        output_plot_path = NULL, num_sig_genes = 20) {


  sdf = JishnuLabTools:::safely_load_obj_from_path(sig_genes)

  er_results = JishnuLabTools:::safely_load_obj_from_path(er_input)

  slide_results = JishnuLabTools:::safely_load_obj_from_path(slide_res)

  x = JishnuLabTools:::safely_load_obj_from_path(xdf)

  y = JishnuLabTools:::safely_load_obj_from_path(ydf)

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

    num_to_plot = ifelse(nrow(temp) < num_sig_genes, nrow(temp), num_sig_genes)
    max_num_to_plot = ifelse(num_to_plot > max_num_to_plot, num_to_plot, max_num_to_plot)

    temp$heights = seq(1, num_to_plot)
    temp$fac = s

    cdf = rbind.data.frame(cdf, temp)
  }

  cdf$loading_anno = ifelse(cdf$A_loading == 1, "*", " ")

  cdf$names_anno = paste0(cdf$names, cdf$loading_anno)

  # get average length of each of the feature names
  avg_feature_length = mean(stringr::str_length(cdf$names))


  plt = cdf %>% ggplot2::ggplot(., aes(x = factor(fac), y = heights, label = names)) +
    ggplot2::geom_text(aes(color = factor(color))) +
    ggplot2::scale_color_manual(values = c("blue", "red"), guide = "none") + theme_void() +
    ggplot2::theme(axis.text.x = element_text(), axis.title.x = element_text(),
          axis.title.y = element_text(angle = 90)) +
    ggplot2::xlab("Latent Factor") + ylab("Genes Associated with Latent Factor Cluster") + ggplot2::ylim(0, num_sig_genes) +
    ggplot2::ggtitle("Genes Associated with Latent Factor")



  if ( !is.null(output_plot_path) ) {


    saveRDS(cdf, paste0(output_plot_path, 'plotSigGenes_data.RDS'))
    ggplot2::ggsave(plot = plt, filename = paste0(output_plot_path, 'plotSigGenes.png'),
                    device = "png", width = 1.5 * length(ks), height = 7)
  }

  return(list("plt" = plt, "plot_df" = cdf))
}
