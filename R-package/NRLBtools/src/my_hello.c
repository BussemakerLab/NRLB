#include <stdio.h>

char* my_hello(void) {
  return("hello!!\n");
}

char* my_string(char *str) {
  return(str);
}

char* my_arg_echo(int argc, char* argv[]) {
  return(argv[0]);
}