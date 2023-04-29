# arglist = list(x_path = list("options" = c("x", "x_path"), "default" = ""),
#                y_path = list("options" = c("y", "y_path"), "default" = ""),
#                out_path = list("options" = c("o", "out_path"), "default" = "/"),
#                run_all = list("options" = c("a", "run_all"), "default" = F))

arg_loader = function(args, default_arglist) {

  arglist = default_arglist
  for(a in arglist) {
    if (length(intersect(a$options, names(args))) > 0) {
      arglist[a] = unlist(args[intersect(a$options, names(args))[1]])
    } else {
      # empty argument; set as default
      arglist[a] = a$default
    }
  }
  return(arglist)
}
