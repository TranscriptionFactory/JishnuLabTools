
get_run_dirs_grep = function(path = ".", filename = "sig_genes.rds") {

  # make filename for searching
  filename_pattern = paste0("(?i)", filename)

  # get everything that ends with r_lambda..._delta... (signature for repeated runs)
  all_runs = grep(x = list.dirs(path),
                  pattern = "lambd(a)?([0-9]+)(\\.)?([0-9]*)_delt(a)?([0-9]+)(\\.)?([0-9]*)$",
                  value = T, ignore.case = T)

  # could also directly just look for ER/SLIDE results
  runs_with_filename = grep(x = list.files(all_runs, full.names = TRUE,recursive = TRUE),
                             pattern = filename_pattern,
                             value = T, ignore.case = T)

  # remove the "sig_genes.rds" so that we know which runs we want to check (we'll need
  # other files from there)
  valid_runs = stringr::str_split_i(runs_with_filename, i = 1,
                                               pattern = filename_pattern)

  # get top level (run) directory for summaries
  run_dir = stringr::str_split_i(valid_runs, i = 1,
                               pattern = "(?<=/)([:alnum:]+_)lambd(a)?([0-9]+)(\\.)?
                               ([0-9]*)_delt(a)?([0-9]+)(\\.)?([0-9]*)")



  return(list(valid_runs = valid_runs, run_directories = run_dir))
}
