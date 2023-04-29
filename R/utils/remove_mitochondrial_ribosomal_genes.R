
remove_mitochondrial_ribosomal_genes = function(df) {

  # check if genes have cluster names e.g C1.RPS prepended
  cluster_names = grepl("^[[:alnum:]].\\d+.HSP|RPS[[:digit:]]|RPL[[:digit:]]|RPLP[[:digit:]]|RPSA|RPS|MT|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|MTND|ATP", colnames(df))

  no_cluster_names = grepl("^HSP|RPS[[:digit:]]|RPL[[:digit:]]|RPLP[[:digit:]]|RPSA|RPS|MT|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|MTND|ATP", colnames(df))

  if (length(cluster_names) > 0) {
    df = df[, -which(cluster_names)]
  } else {
    df = df[, -which(no_cluster_names)]
  }
  return(df)
}
