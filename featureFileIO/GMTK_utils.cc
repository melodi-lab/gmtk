
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <ctype.h>
#include "GMTK_utils.h"


gmtk_off_t
fileSize(FILE *f) {

  gmtk_off_t s;

  if (gmtk_fseek(f,(gmtk_off_t)0,SEEK_END) != 0)
    error("fileSize: Can't skip to end of file\n");
  s = gmtk_ftell(f);

  rewind(f);
  return s;
}

char *
skip_whitespace(char *s) {

  char *tmp,*endp;
  tmp = s;
  endp = &(s[strlen(s)-1]);

  while (tmp != endp) {
    if (!(isspace(*tmp)))
      return tmp;
  }
  return tmp;
}

char *
skip_chars(char *s) {


  char *endp = &s[strlen(s)];
  char *tmp = s;

  while (tmp != endp) {
    if (isspace(*tmp))
      return tmp;
  }
  return tmp;
}

int
emptyline(char *p1, char *p2) {

  char *p = p1;

  if (p1 == p2)
    return 0;

  while (p != p2 ) {
    if (!(isspace(*p))) {
      return 0;
    }
    p++;
  }
  return 1;
}





