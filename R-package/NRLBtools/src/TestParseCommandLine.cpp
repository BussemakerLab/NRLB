#include <Rcpp.h>
using namespace Rcpp;

#include <string.h>

extern "C" {
  #include "my_hello.h"
}

//' Run 'LogoGenerator' tool from REDUCE Suite
//' 
//' @param argv A commandline option string
//' @return An optional error code
//' @export
// [[Rcpp::export]]
void TestParseCommandLine(std::string commandline) {
  
  char toolname[100] = "LogoGenerator";
    
  int argc = 0;
  char *argv[100];

  // argv[0] should contain tool name
  argv[argc++] = toolname;
  
  // split commandline by whitespace to further populate argv[]
  char str[1000];
  char *token;
  char *rest;
  std::strcpy(str, commandline.c_str());
  rest = str;
  while ((token = strtok_r(rest, " ", &rest))) {
    argv[argc++] = token;
    Rprintf("token='%s' str='%s' rest='%s'\n", token, str, rest);
  }
  Rprintf(my_arg_echo(argc, argv));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*** R
TestParseCommandLine("--arg1 --arg2 -arg3==val1")
*/
