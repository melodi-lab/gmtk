#source this file to add all the binaries and scripts to your path,
#and the required modules to the module paths of perl and python.
#
# Arthur 1/1/09

PREFIX=/cworkspace/ifp-32-1/hasegawa/programs/
export PATH=$PREFIX/gmtk/bin.`uname`:$PATH
export PATH=$PREFIX/htk:$PATH
export PATH=$PREFIX/sctk-2.3/bin.`arch`:$PATH
export PATH=$PREFIX/srilm/bin.`uname`:$PATH
export PATH=$PREFIX/bin.`arch`:$PATH
#export PATH=$PREFIX/bin.`uname`:$PATH
export PATH=$PREFIX/quicknet/bin.`arch`:$PATH


SCRPREFIX=$PREFIX/scripts
export PATH=$SCRPREFIX/parallel:$SCRPREFIX/gmtk:$PATH
export PERL5LIB=$PERL5LIB:$SCRPREFIX/perlLib:$SCRPREFIX/perlLibThirdParty
export PYTHONPATH=$SCRPREFIX/pythonLib
