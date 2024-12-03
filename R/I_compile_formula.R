
.compile_formula <- function(form){
  env <- environment(form)
  dec <- terms(form, keep.order = TRUE)
  ter <- attr(dec, "term.labels")

  if(length(ter)>1){stop("Too many effects. Only one is permitted")}
  if(length(which(startsWith(ter, "s_n")))==0){stop("There is no nested effect in the formula")}

  eval(parse(text = ter), envir = env)
}
