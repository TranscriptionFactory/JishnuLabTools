library(tidyverse)
library(qgraph)

# here is a simple function
create_corr_network = function(x_data, slide_res) {
  # slide res should be the output from SLIDE + SLIDEHelper that has the fields:
  # slide_res$SLIDE_res$marginal_vars - vector with marginal variables
  # slide_res$feature_res - list of DFs with top genes from each LF

  x = cor(x)

}


create_corr_network = function(x, y, out_dir, plot_name, gene_names) {

  # x - sample x feature matrix
  # y - labels (assuming binary)
  # out_dir - where you want the plot saved
  # plot_name - title for plot
  # gene_names - vector with sig genes you want to select
  setwd(out_dir)

  x_gene = as.matrix(x[, gene_names])

  # uncomment if you want the nodes colored by AUC
  # col_auc = round(apply(x_gene, 2, function(xs) glmnet:::auc(as.matrix(y), as.matrix(xs))), 2)
  # temp_cols = ifelse(col_auc > 0.55, "#40006D", ifelse(col_auc < 0.45, "#59A14F", "lightgray"))

  x_temp = cor(x_gene)

  pl = qgraph(x_temp, filename= plot_name,
              layout = "spring", threshold=0.4, repulsion=0.1,
              labels = colnames(x_temp),
              # color = temp_cols, # uncomment if you want the nodes colored by AUC
              title = plot_name,
              label.scale.equal=FALSE,label.prop=0.95,shape="ellipse",
              posCol="salmon", negCol="skyblue",filetype='pdf',height=5,width=7)
}
