#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
  #include "LogoGenerator_Rinterface.h"
}

// [[Rcpp::export]]
void LogoGenerator_Cpp(std::string output, 
                       std::string file, 
                       std::string logo, 
                       std::string title, 
                       double ymin, 
                       double ymax) {
  char output_cstr[100], file_cstr[100], logo_cstr[100], title_cstr[100];
  std::strcpy(output_cstr, output.c_str());
  std::strcpy(file_cstr, file.c_str());
  std::strcpy(logo_cstr, logo.c_str());
  std::strcpy(title_cstr, title.c_str());
  
  LogoGenerator_C(output_cstr, 
                  file_cstr, 
                  logo_cstr, 
                  title_cstr, 
                  ymin, 
                  ymax);
}

/*** R
LogoGenerator(output="output", file="file", log="logo", 
              title="title", ymin=0.0, ymax=1.0)
*/