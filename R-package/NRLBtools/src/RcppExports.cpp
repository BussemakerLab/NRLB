// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// LogoGenerator_Cpp
void LogoGenerator_Cpp(std::string output, std::string file, std::string logo, std::string title, double ymin, double ymax);
RcppExport SEXP _NRLBtools_LogoGenerator_Cpp(SEXP outputSEXP, SEXP fileSEXP, SEXP logoSEXP, SEXP titleSEXP, SEXP yminSEXP, SEXP ymaxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type output(outputSEXP);
    Rcpp::traits::input_parameter< std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type logo(logoSEXP);
    Rcpp::traits::input_parameter< std::string >::type title(titleSEXP);
    Rcpp::traits::input_parameter< double >::type ymin(yminSEXP);
    Rcpp::traits::input_parameter< double >::type ymax(ymaxSEXP);
    LogoGenerator_Cpp(output, file, logo, title, ymin, ymax);
    return R_NilValue;
END_RCPP
}
// TestParseCommandLine
void TestParseCommandLine(std::string commandline);
RcppExport SEXP _NRLBtools_TestParseCommandLine(SEXP commandlineSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type commandline(commandlineSEXP);
    TestParseCommandLine(commandline);
    return R_NilValue;
END_RCPP
}
// TestParseStringVectorArgs
void TestParseStringVectorArgs(Rcpp::StringVector arguments);
RcppExport SEXP _NRLBtools_TestParseStringVectorArgs(SEXP argumentsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type arguments(argumentsSEXP);
    TestParseStringVectorArgs(arguments);
    return R_NilValue;
END_RCPP
}
// timesTwo
NumericVector timesTwo(NumericVector x);
RcppExport SEXP _NRLBtools_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NRLBtools_LogoGenerator_Cpp", (DL_FUNC) &_NRLBtools_LogoGenerator_Cpp, 6},
    {"_NRLBtools_TestParseCommandLine", (DL_FUNC) &_NRLBtools_TestParseCommandLine, 1},
    {"_NRLBtools_TestParseStringVectorArgs", (DL_FUNC) &_NRLBtools_TestParseStringVectorArgs, 1},
    {"_NRLBtools_timesTwo", (DL_FUNC) &_NRLBtools_timesTwo, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_NRLBtools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
