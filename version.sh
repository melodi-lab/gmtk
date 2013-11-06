#! /bin/sh

if test -d .hg ; then
  echo GMTK release: `cat RELEASE`
  echo Mercurial revision id: `hg id`
  echo Mercurial check-in date: `hg parents --template '{date|date}'`
  echo Mercurial branch: `hg branch`
else
  cat version.info
fi
