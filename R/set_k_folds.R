#' @export
set_k_folds = function(er_input){

  # load the y from input
  y = JishnuLabTools:::safely_load_obj_from_path(er_input$y_path)
  #as.matrix(read.csv(er_input$y_path), row.names = 1)
  # sample size
  n = min(length(y), nrow(y))

  # default is 10 folds unless we have small sample size
  k = 10
  if (n < 15) {
    # LOOCV
    k = n
  } else if (n < 30) {
    k = 5
  }
  return(k)
}
