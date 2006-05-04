plotMD <- function(object, type = c("coverage", "density"), ...) {
  if (type == "coverage")
    plot.eiMDcoverage(object, ...) 
  else if (type == "density")
    plot.eiMDdensity(object, ...)
  else
   stop(paste(type, "plot for class 'eiMD' not supported."))
}
