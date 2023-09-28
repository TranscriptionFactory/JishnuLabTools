## this code makes correlation networks ##
# load packages
library(tidyverse)
library(qgraph) ## for making the network

create_corr_network = function(er_run_path, x_path, y_path) {

  x = read.csv(x_path, row.names = 1)

  y = read.csv(y_path)

  plot_sig_genes_path = list.files(er_run_path, recursive = T,
                                   pattern = "plotSigGenes_data.RDS",
                                   full.names = T)

  sig_genes = readRDS(plot_sig_genes_path)

  out_dir = paste0(er_run_path, "/LF_corr_network")
  if ( !dir.exists(out_dir) ) {
    dir.create(out_dir)
  }

  for (f in unique(sig_genes$fac)) {

        x_gene = x[, (sig_genes %>% filter(fac == f))$names]

    col_auc = round(apply(x_gene, 2, function(x) glmnet:::auc(as.numeric(y), as.matrix(x))), 2)

    temp_cols = ifelse(col_auc > 0.55, "#40006D", ifelse(col_auc < 0.45, "#59A14F", "lightgray"))

    x_temp = cor(x_gene)

    # change wd cause weird
    setwd(out_dir)

    pl = qgraph(x_temp, filename= paste0(f, "_corr_network"),
                layout = "spring", threshold=0.5, repulsion=0.1,
                labels = colnames(x_temp),
                # title = paste0(comp, " ", color_code),
                label.scale.equal=FALSE,label.prop=0.95,shape="ellipse",
                posCol="salmon", negCol="skyblue",filetype='pdf',height=3,width=3)

  }
}
