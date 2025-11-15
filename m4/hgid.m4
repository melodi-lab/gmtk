AC_DEFUN([AC_HG_ID],
[ AC_CHECK_FILE([$srcdir/../.git/config],[
  AC_CHECK_PROG([HAVE_HG],[git],[yes])
  if test x$HAVE_HG = xyes; then
    HGREVID=`git rev-parse HEAD`
    HGDATE=`git log --no-walk --first-parent --format='format:%aD' HEAD`
    printf "#ifndef HGID\n#define HGID \"$HGREVID\"\n#define HGDATE \"$HGDATE\"\n#endif\n" > hgstamp.h
    AC_DEFINE([HAVE_HG_H],[1],[Live Mercurial revision id available?])
  fi
],[
  AC_CHECK_FILE([hgstamp.h],[
    AC_DEFINE([HAVE_HG_H],[1],[Distributed Mercurial revision id header available?])
  ],[
    AC_DEFINE([HGID],"UNKNOWN",[Unknown Mercurial revision id - broken distribution?])
    AC_DEFINE([HGDATE],"UNKOWN",[Unkown Mercurial checkin date - broken distribution?])
  ])
])
AM_CONDITIONAL([HAVE_HG],[test  x$HAVE_HG = xyes])
])
