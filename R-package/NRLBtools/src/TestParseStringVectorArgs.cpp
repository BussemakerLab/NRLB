#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
  #include "my_hello.h"
}

//' Run 'LogoGenerator' tool from REDUCE Suite
//' 
//' @param argv A commandline option string
//' @return An optional error code
//' @export
// [[Rcpp::export]]
void TestParseStringVectorArgs(Rcpp::StringVector arguments) {
  
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
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*** R
 TestParseStringVectorArgs(c("arg1", "arg2", "arg3=val3"))
*/
