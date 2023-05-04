arglist = list(x_path = list(options = c("x", "x_path"), default = ""),
               y_path = list(options = c("y", "y_path"), default = ""),
               data_path = list(options = c("d", "data", "data_path"),
                                default = ""),
               yaml_path = list(options = c("yaml", "yaml_path", "p"),
                                default = ""),
               out_path = list(options = c("o", "out_path"), default = "/"),
               coarse_grid = list(options = c("c", "coarse_grid", "coarseGrid"),
                                  default = F))

arg_loader = function(args, default_arglist) {

  arglist = default_arglist
  for(a in arglist) {
    if (any(sapply(a$options, function(x) any(stringr::str_detect(names(args),
                                                              pattern = x))))) {
      #length(intersect(a$options, names(args))) > 0) {
      arglist[[a]] = unlist(args[intersect(a$options, names(args))[1]])
    } else {
      # empty argument; set as default
      arglist[[a]] = a$default
    }
  }
  return(arglist)
}
