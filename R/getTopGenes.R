# https://github.com/Hanxi-002/SLIDEHelper/blob/main/R/GetTopFeatures.R
# #' Output the top features in each latent factors
# #'
# #' Output the top features in each latent factors by taking a union of features with highest A loadings and spearman correlations or univariate AUCs
# #'
# #' @param x_path a string that points to the x matrix.
# #' @param y_path a string that points to the y vector.
# #' @param er_path a string that points to the final er result as an RDS object.
# #' @param out_path a string that points to a folder where txt files of each latent factors will be outputted.
# #' @param SLIDE_res the object outputted by the runSLIDE function. 
# #' @param num_top_feats the number of top features to get from each criteria such as A loadings or AUCs.
# #' @param condition either "corr" or "auc", indicates which one to use to get top features. 
# #' @return The SLIDE_res object with feature results added to one of the slots.
# #' @export

# # #' @export
# GetTopFeatures <- function(x, y, er_res, SLIDE_res, num_top_feats = 20){
#   x <- as.matrix(utils::read.csv(er_res$x_path, row.names = 1))
#   y <- as.matrix(utils::read.csv(er_res$y_path, row.names = 1))

#   condition = er_res$eval
  
#   ks <- union(union(unique(SLIDE_res$interaction$p1), unique(SLIDE_res$interaction$p2)), unique(SLIDE_res$marginal_vals))
#   if (is.null(ks) == TRUE){stop('The SLIDE_res input is not formatted correctly. Please re-run the runSLIDE function...')}
#   A <- er_res$A[, ks]
#   gene_names <- colnames(x)
  
#   temp <- NULL
#   for (i in 1:ncol(A)){
#     AUCs <- c()
#     signs <- c()
#     corrs <- c()
#     idx <- which(A[, i] != 0)
#     A_loading <- abs(A[, i][-which(A[, i] == 0)])
#     names <- gene_names[idx]
    
#     for (j in 1:length(idx)) { ## loop through variables with big enough loadings
#       corr <- cor(x[,idx[j]],y,method = "spearman")
#       sign <- sign(cor(x[,idx[j]],y,method = "spearman"))
      
#       if (condition == "auc"){
#         AUC <- auc(y, x[, idx[j]])
#         AUCs <- c(AUCs, AUC)
#       }
      
#       corrs <- c(corrs, corr)
#       signs <- c(signs, sign)
#     }
#     color <- recode(signs, "-1" = "Blue", "1"= "Red")
#     if (condition == "auc"){
      
#       df <- data.frame(names, A_loading, AUCs, corrs, color)
#       df <- df[order(abs-df$A_loading), ]
#       top <- df[1:num_top_feats, ]
      
#       df <- df[order(-df$AUCs), ]
#       bot <- df[1:num_top_feats, ]
#     }else if(condition == "corr"){
      
#       df <- data.frame(names, A_loading, corrs, color)
#       df <- df[order(-df$A_loading), ]
#       top <- df[1:num_top_feats, ]
      
#       df <- df[order(-abs(df$corrs)), ]
#       bot <- df[1:num_top_feats, ]
#       final <- unique(rbind(top, bot))
#     }
    
#     temp[[i]] <- final
#     write.table(final, file = paste0(out_path, "/gene_list_",
#                                      colnames(A)[i], ".txt"), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
#   }
#   names(temp) <- colnames(A)
#   SLIDE_res$feature_res <- temp
#   return(SLIDE_res)
# }
getTopGenes = function(er_res, ks, x, y) {
  A <- as.matrix(er_res$A[, ks])
  gene_names <- colnames(x)
  temp <- NULL
  for (i in 1:ncol(A)){
    AUCs <- c()
    signs <- c()
    corrs <- c()
    idx <- which(A[, i] != 0)
    A_loading <- abs(A[, i][-which(A[, i] == 0)])
    names <- gene_names[idx]
    for (j in 1:length(idx)) { ## loop through variables with big enough loadings
      corr <- cor(x[,idx[j]],y,method = "spearman")
      sign <- sign(cor(x[,idx[j]],y,method = "spearman"))
      AUC <- glmnet:::auc(y, x[, idx[j]])
      AUCs <- c(AUCs, AUC)
      corrs <- c(corrs, corr)
      signs <- c(signs, sign)
    }
    color <- recode(signs, "-1" = "Blue", "1"= "Red")
    df <- data.frame(names, A_loading, AUCs, corrs, color)

    num_genes = 10
#     if( nrow(df) < 10 ) {
#       num_genes = nrow(df) / 2
#     }
    df <- df[order(-df$A_loading), ]
    top <- df[1:num_genes, ]
    df <- df[order(-df$AUCs), ]
    bot <- df[1:num_genes, ]
    final <- rbind(top, bot)
    temp[[i]] <- final
  }
  return(temp)
}
