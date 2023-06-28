
plotSigGenes_SLIDEHelper = function(slide_results, output_plot_path = NULL,
                        annotate_anchors = F) {

  slide_marginals = slide_results$feature_res

  # this is for the upper y limit in the plot
  max_num_genes_in_any_lf = 0

  sg_df = data.frame()

  # make results into plottable dataframe
  for (lf in 1:length(slide_marginals)) {
    lf_df = list2DF(slide_marginals[[lf]])
    lf_df$lf_num = as.numeric(stringr::str_replace(names(slide_marginals)[lf],
                                    pattern = "Z", replacement = ""))
    sg_df = rbind(sg_df, lf_df)

    max_num_genes_in_any_lf = ifelse(nrow(lf_df) > max_num_genes_in_any_lf,
                                     nrow(lf_df),
                                     max_num_genes_in_any_lf)
  }

  sg_plot_df = data.frame()
  for (lf in unique(sg_df$lf_num)) {
    # get genes for this latent factor
    lf_temp = sg_df[sg_df$lf_num == lf, ]

    # add the heights to each gene
    lf_temp$plot_height = seq(1, nrow(lf_temp))

    sg_plot_df = rbind(sg_plot_df, lf_temp)
  }

  # if want to have anchors annotated with *
  if (annotate_anchors) {
    sg_plot_df$loading_anno = ifelse(sg_plot_df$A_loading == 1, "*", " ")
    sg_plot_df$names = paste0(sg_plot_df$names, sg_plot_df$loading_anno)
  }


  plt = sg_plot_df %>% ggplot2::ggplot(., aes(x = factor(lf_num), y = plot_height, label = names)) +
    ggplot2::geom_text(aes(color = factor(color))) +
    ggplot2::scale_color_manual(values = c("blue", "red"), guide = "none") + theme_void() +
    ggplot2::theme(axis.text.x = element_text(), axis.title.x = element_text(),
                   axis.title.y = element_text(angle = 90)) +
    ggplot2::xlab("Significant Latent Factor") +
    ylab("Genes Associated with Significant Latent Factors") +
    ggplot2::ylim(0, max_num_genes_in_any_lf) +
    ggplot2::ggtitle("Genes Associated with Significant Latent Factors")


  if ( !is.null(output_plot_path) ) {

    saveRDS(sg_plot_df, paste0(output_plot_path, 'plotSigGenes_data.RDS'))
    ggplot2::ggsave(plot = plt, filename = paste0(output_plot_path, 'plotSigGenes.png'),
                    device = "png",
                    width = 1.5 * length(unique(sg_plot_df$lf_num)), height = 7)
  }

  return(list("plt" = plt, "plot_df" = sg_plot_df))
}
