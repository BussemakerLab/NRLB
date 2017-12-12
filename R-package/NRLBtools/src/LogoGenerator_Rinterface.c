#include <stdio.h>
#include <string.h>
#include "REDUCE_Suite.h"


void LogoGenerator_C(char *output, 
                     char *file, 
                     char *logo, 
                     char *title, 
                     double ymin, 
                     double ymax) {

  args_logo_generator args;

  set_my_globals("LogoGenerator");
  set_logo_defaults(&args);
  
  strcpy(args.outdir, output);
  strcpy(args.file,   file);
  strcpy(args.logo,   logo);
  strcpy(args.title,  title);
  strcpy(args.style,  "ddG");
  strcpy(args.format, "PNG");
  
  printf("args.outdir='%s'\n", args.outdir);
  Gvars.PRGLOG = stdout;
  printf("Gvars.PRGLOG='%p'\n", Gvars.PRGLOG);

  logo_psam(&args);
}

int main(int argc, char **argv) {
  
  LogoGenerator_C("/Users/HJB/GitHub/NRLB/R-package/NRLBtools/src/junkdir",
                  "/Users/HJB/REDUCE-Suite-v2.2/examples/MatrixREDUCE/spellman-alpha/results/psam_001.xml",
                  "logofile",
                  "this is the title", 
                  0.0, 1.0);
  return(0);
  
}