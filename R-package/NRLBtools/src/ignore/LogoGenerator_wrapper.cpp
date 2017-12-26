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
void LogoGenerator_alt1() {
  std::string message = "blahhh!";
  Rprintf(message);
}

//' Run 'LogoGenerator' tool from REDUCE Suite
//' 
//' @param argv A commandline option string
//' @return An optional error code
//' @export
// [[Rcpp::export]]
void LogoGenerator(std::string commandline) {
  
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

//' Run 'LogoGenerator' tool from REDUCE Suite
//' 
//' @param argv A commandline option string
//' @return An optional error code
//' @export
void LogoGenerator_alt3(Rcpp::StringVector arguments) {
  
  int argc = arguments.size();
  Rprintf("argc=%d\n", argc);
  
  //std::vector<std::string> argv2(argc);
  char argv[1000][1000];
  int index;
  std::string str;
  for (index = 0; index < argc; index++){
    str = arguments(index);
    std::strcpy(argv[index], str.c_str());
    Rprintf("argv[%d]='%s'\n", index, argv[index]);
    
    //argv2[index] = Rcpp::as< std::string > (arguments(index));
    //Rprintf("argv2[%d]='%s'\n", index, argv2[index].c_str());
  }

  Rprintf(my_hello());
  Rprintf(my_string(argv[1]));

  // this version not working; problem passing vector of strings as char *argv[] to C-level main() function
  
  //char** test = argv;
  //Rprintf(my_arg_echo(argc, &argv[0]));
  
  //C_main_LogoGenerator(argc, argv);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*
/*** R
LogoGenerator("bla")
*/
