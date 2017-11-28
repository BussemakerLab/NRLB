#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
    #include "my_hello.h"
    //#include "LogoGenerator_main.h"
    //int main(int argc, char *argv[]);
}

//' Run 'LogoGenerator' tool from REDUCE Suite
//' 
//' @param argv A commandline option string
//' @return An optional error code
//' @export
// [[Rcpp::export]]
void LogoGenerator(StringVector arguments) {
  
  int argc = arguments.size();
  Rprintf("argc=%d\n", argc);
  
  char argv[1000][1000];
  int index;
  std::string str;
  for (index = 0; index < argc; index++){
    str = arguments(index);
    strcpy(argv[index], str.c_str());
    Rprintf("argv[%d]='%s'\n", index, argv[index]);
  }
  //strcpy(argv[argc], "");
  Rprintf(my_hello());
  
  //main(argc, (char **) argv);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
LogoGenerator(c("LogoGenerator", "--file=blahhhhh", "--verbose"))
*/
