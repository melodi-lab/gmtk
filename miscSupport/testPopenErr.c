
#include <stdlib.h>
#include <stdio.h>
#include "popen_err.h"

#define MAX_COMMAND_LEN 1024

int
main(int argc, char *argv[]) {
  FILE *f = popen_err("./a.out blah foo bar glort", "r", "prefix: ");
  char s[MAX_COMMAND_LEN];
  if (!f) {
    perror("open_err");
    exit(EXIT_FAILURE);
  }
  while (1) {
    fgets(s, MAX_COMMAND_LEN-1, f);
    if (feof(f) || ferror(f)) break;
    printf("%s", s);
  } 
  if (pclose_err(f) == -1) {
    perror("pclose_err");
    exit(EXIT_FAILURE);
  }
  exit(0);
}
