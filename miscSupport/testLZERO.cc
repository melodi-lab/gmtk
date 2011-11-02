
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "logp.h"

int
main(int argc, char *argv[]) {
  logpr a(0.0);
  logpr b(0.0);
  logpr c = a + b;
  logpr d = a - b;
  if (c.zero() && d.zero())
    return 0;
  else
    return 5;
}
