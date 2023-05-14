
# runs: vector of file paths
# filename: common name shared among all output files that you want (e.g all
# ER pipeline3 outputs begin with 'final')
load_output_from_runs = function(runs, filename) {
                                #function_to_read_files = JishnuLabTools:::safely_load_obj_from_path) {

  runs_with_desired_file = lapply(runs, function(x)
    list.files(x, full.names = T)[stringr::str_which(list.files(x),
                                                     pattern = filename)])

  # remove any entries that didn't find a match
  for (l in 1:length(runs_with_desired_file)) {
    if (length(runs_with_desired_file[[l]]) == 0) {
      runs_with_desired_file[[l]] = NULL
    } else if (length(runs_with_desired_file[[l]]) > 1) {
      # if theres more than one match, just pick first
      runs_with_desired_file[[l]] = runs_with_desired_file[[l]][1]
    }
  }

  if ( length(runs_with_desired_file) == 0 ) {
    return(NULL)
  } else {
    # read files using function passed
    res = lapply(runs_with_desired_file, function(x) JishnuLabTools:::safely_load_obj_from_path(x))#function_to_read_files(x))
    return( res )
  }
}
