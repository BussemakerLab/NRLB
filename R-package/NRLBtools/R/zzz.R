.onLoad <- function(libname, pkgname){
  packageStartupMessage("Welcome to the NRLBtools package")
  Sys.setenv(REDUCE_SUITE=system.file(package="NRLBtools"))
}