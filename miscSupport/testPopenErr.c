
#include <stdlib.h>
#include <stdio.h>
#include "popen_err.h"

int
main(int argc, char *argv) {
  FILE *f = popen_err("./a.out", "r", "prefix: ");
  if (!f) {
    perror("open_err");
    exit(EXIT_FAILURE);
  }
  char s[1024];
  for (; !feof(f); ) {
    fgets(s, 1023, f);
    printf("%s", s);
  }
  if (pclose_err(f) == -1) {
    perror("pclose_err");
    exit(EXIT_FAILURE);
  }
  exit(0);
}
