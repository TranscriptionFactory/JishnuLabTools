## this code makes correlation networks ##
# load packages
library(tidyverse)
library(qgraph) ## for making the network




create_corr_network = function(path, plot_label = NULL) {

  yaml_input = yaml::yaml.load_file(path)

  dir_name = yaml_input$out_path

  x_path = yaml_input$x_path
  y_path = yaml_input$y_path

  x_mat = as.matrix(read.csv(x_path, row.names = 1, check.names = F))
  y_mat = as.matrix(read.csv(y_path, row.names = 1, check.names = F))


############### uncomment if you want to use the ER path instead
  check_for_file = function(path, file_pattern) {

    f = list.files(path, recursive = T, full.names = T, pattern = file_pattern)
    if( length(f) > 0 ) {
      return(f[1])
    } else {
      cat("\n Folder must have exactly one file with desired name: ", pattern, "\n")
      return(NULL)
    }
  }

  get_x_and_y = function(path) {
    x = check_for_file(path, file_pattern = "x.csv")
    y = check_for_file(path, file_pattern = "y.csv")

    load_matrix_csv = function(path) {
      return(read.csv(path, row.names = 1))
    }

    if (all(!is.null(x), !is.null(y))) {
      return(cbind.data.frame(load_matrix_csv(y), load_matrix_csv(x)))
    } else {
      df = check_for_file(path, file_pattern = df_filename_RDS_pattern)
      if (!is.null(df)) {
        return(readRDS(df))
      } else {
        return(NULL)
      }
    }
  }


  df = get_x_and_y(path)
  if (is.null(df)) {
    cat("\n Failed to load dataframes. Check path \n")
    return()
  } else {
    names(df)[1] = "y"
  }

  sig_genes_data = check_for_file(dir_name, file_pattern = "plotSigGenes_data.RDS")
  sig_genes_data = readRDS(sig_genes_data)

  all_LFs = c()

  if ( !is.null(sig_genes_data) ) {
    all_LFs = stringr::str_to_upper(sig_genes_data$names)

  } else {
    cat("\n Couldn't find genes in significant latent factors. Check path. \n")
    return()
  }

  color_code = stringr::str_to_lower(sig_genes_data$color)

  x = df[, -1]

  if ( dim(x)[2] <= 1) {
    cat("\n Latent factor genes not found in X data. Check path \n")
    return()
  }

    dir_name = paste0(dir_name, "/LF_correlation_plots/")

  if ( !dir.exists(dir_name) ) {
    dir.create(dir_name)
  }

  original_wd = getwd()

  for ( LF in unique(sig_genes_data$fac) ) {
    sg_temp = ( sig_genes_data %>% filter(fac == LF) )

    x_cor = x[, sg_temp$names]

    x_temp = cor(as.matrix(x_cor))

    col_auc = round(apply(x_gene, 2, function(xs) glmnet:::auc(as.matrix(y), as.matrix(xs))), 2)
    
    node_color = ifelse(col_auc > 0.55, "#40006D", ifelse(col_auc < 0.45, "#59A14F", "lightgray"))
    
    node_shape = ifelse(is.na(sg_temp$color), "circle", ifelse(sg_temp$color == "Red", "triangle", "square"))

    # change wd cause weird
    setwd(dir_name)
    plot_title = ifelse(is.null(plot_label), "", plot_label)

    # plot_title = paste0("LF #", LF, "\n", plot_title, "\nTriangle = Up\nSquare = Down")
    plot_title = paste0("LF #", LF)

    # this won't get plotted. the only important part is to have the 'groups' argument be whatever you want to group
    # by - e.g if you want the shapes clustered together or if you want nodes with a specific color clustered together
    pl_shape = qgraph(x_temp, filename = LF, groups = node_shape,
                legend.mode = "groups",
                color = clust_colors,
                shape = node_shape,
                layout = "spring",
                # minimum=0.1, # use minimum instead of threshold
                labels = colnames(x_temp),
                label.scale.equal=TRUE,label.prop=0.8,shape="ellipse",
                DoNotPlot = TRUE,
                posCol="#51087E", negCol="#59A14F",filetype='pdf',height=10,width=10)

    # change the 'groups' argument here to whatever you want the legend to show
    pl = qgraph(x_temp, filename = paste0("celltype_LF_shape_cluster_fade_", LF), groups = clust_names_list,
                color = color_names_vec,
                shape = node_shape,
                edge.width = 1,
                node.width = 0.75,
                node.height = 0.75, fade = TRUE, # turn off fade if you want the edge color to be constant
                minimum=0.25,
                labels = colnames(x_temp),
                # title = plot_title,
                layout = pl_shape$layout, # use the plot layout above
                label.scale.equal=TRUE,label.prop=0.95,shape="ellipse",
                posCol="#51087E", negCol="#59A14F",filetype='pdf',height=5,width=8)
  }

  setwd(original_wd)
}


