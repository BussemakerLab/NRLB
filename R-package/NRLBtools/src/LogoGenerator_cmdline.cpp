#include <Rcpp.h>
using namespace Rcpp;

#include <string.h>

extern "C" {
  #include "LogoGenerator_Rinterface.h"
}

// [[Rcpp::export(.LogoGenerator)]]
void LogoGenerator_cmdline(std::string argstr) {
  
  char toolname[100] = "LogoGenerator";
    
  int argc = 0;
  char *argv[100];

  // argv[0] should contain tool name
  argv[argc++] = toolname;
  
  // split commandline by whitespace to further populate argv[]
  char str[1000];
  char *token;
  char *rest;
  std::strcpy(str, argstr.c_str());
  rest = str;
  while ((token = strtok_r(rest, " ", &rest))) {
    argv[argc++] = token;
  }
  
  // run LogoGenerator
  LogoGenerator_main(argc, argv);
}