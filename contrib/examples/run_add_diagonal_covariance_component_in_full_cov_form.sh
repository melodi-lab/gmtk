#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#
exe_name=$0
exe_dir='add_diagonal_covariance_component_in_full_cov_form_buildfiles'
MCRROOT="/usr/nikola/pkgs/matlab"
MWE_ARCH="glnxa64" ;
if [ "$MWE_ARCH" = "sol64" ] ; then
	LD_LIBRARY_PATH=.:/usr/lib/lwp:${MCRROOT}/runtime/glnxa64 ; 
else
	LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
fi
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
if [ "$MWE_ARCH" = "maci" -o "$MWE_ARCH" = "maci64" ]; then
	DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/System/Library/Frameworks/JavaVM.framework/JavaVM:/System/Library/Frameworks/JavaVM.framework/Libraries;
else
	MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
fi
XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
export LD_LIBRARY_PATH;
export XAPPLRESDIR;
${exe_dir}/add_diagonal_covariance_component_in_full_cov_form $*
exit
