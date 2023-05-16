#' @export
remove_mitochondrial_ribosomal_genes = function(df) {

  # check if genes have cluster names e.g C1.RPS prepended
  cluster_names = which(grepl("^[[:alnum:]].\\d+.RPL|RPS[[:digit:]]|RPL[[:digit:]]|RPLP[[:digit:]]|RPSA|RPS|MT|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|MTND|ATP", stringr::str_to_upper(colnames(df))) == T)

  no_cluster_names = which(grepl("^RPL|RPS[[:digit:]]|RPL[[:digit:]]|RPLP[[:digit:]]|RPSA|RPS|MT|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|MTND|ATP", stringr::str_to_upper(colnames(df))) == T)

  if (length(cluster_names) == 0 && length(no_cluster_names) == 0) {
    return(df)
  } else if (length(cluster_names) > 0) {
    df = df[, -cluster_names]
  } else {
    df = df[, -no_cluster_names]
  }
  return(df)
}
