x_temp = cor(x_mat)

# change wd to where you want to saveuse weird
setwd(out_dir)

# will save as out_dir / corr_network.pdf
pl = qgraph(x_temp, filename= "corr_network",
            layout = "spring", threshold=0.5, repulsion=0.1,
            labels = colnames(x_temp),
            # title = paste0(comp, " ", color_code),
            label.scale.equal=FALSE,label.prop=0.95,shape="ellipse",
            posCol="salmon", negCol="skyblue",filetype='pdf',height=3,width=3)
