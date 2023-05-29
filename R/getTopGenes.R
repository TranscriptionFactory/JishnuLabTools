
#' @export
getTopGenes = function(er_res, ks, x, y, num_genes = 20) {
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
    color <- dplyr::recode(signs, "-1" = "Blue", "1"= "Red")
    df <- data.frame(names, A_loading, AUCs, corrs, color)

    num_genes = ifelse(round(num_genes/2) > 1, round(num_genes/2), 10)

    df <- df[order(-df$A_loading), ]
    top <- df[1:num_genes, ]
    df <- df[order(-df$AUCs), ]
    bot <- df[1:num_genes, ]
    final <- unique(rbind(top, bot))
    temp[[i]] <- final
  }
  return(temp)
}
