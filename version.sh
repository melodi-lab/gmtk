#! /bin/sh

if test -d .git ; then
  echo GMTK release: `cat RELEASE`
  echo git revision id: `git rev-parse HEAD`
  echo git check-in date: `git log --no-walk --first-parent --format='format:%aD' HEAD`
  echo git branch: `git branch`
else
  cat version.info
fi
