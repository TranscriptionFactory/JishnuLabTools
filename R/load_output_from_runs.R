
# runs: vector of file paths
# filename: common name shared among all output files that you want (e.g all
# ER pipeline3 outputs begin with 'final')
load_output_from_runs = function(runs, filename) {
                                #function_to_read_files = JishnuLabTools:::safely_load_obj_from_path) {

  res_paths = get_run_dirs_grep

  if ( length(res_paths$valid_runs) == 0 ) {
    return(NULL)
  } else {
    # read files using function passed
    res = lapply(res_paths$valid_runs, function(x) JishnuLabTools:::safely_load_obj_from_path(x))
    return( res )
  }
}
