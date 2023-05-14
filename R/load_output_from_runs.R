
# runs: vector of file paths
# filename: common name shared among all output files that you want (e.g all
# ER pipeline3 outputs begin with 'final')
load_output_from_runs = function(runs, filename,
                                function_to_read_files = function(x) JishnuLabTools:::safely_load_obj_from_path(x)) {

  runs_with_desired_file = sapply(runs, function(x)
    list.files(x, full.names = T)[stringr::str_which(list.files(x),
                                                     pattern = filename)])
  if ( length(runs_with_desired_file) == 0 ) {
    return(FALSE)
  } else {
    # read files using function passed
    res = lapply(runs_with_desired_file, function(x) function_to_read_files(x))
    return( res )
  }
}
