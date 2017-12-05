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

  strcpy(args.outdir, output);
  strcpy(args.file,   file);
  strcpy(args.logo,   logo);
  strcpy(args.title,  title);
  strcpy(args.style,  "ddG");
  strcpy(args.format, "PNG");
  
  printf("args.outdir='%s'\n", args.outdir);
}