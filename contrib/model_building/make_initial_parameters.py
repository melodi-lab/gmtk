#! /usr/bin/env /g/ssli/transitory/sking/python/Python-2.3.3/python

# ---------------------------------------------------------------
#
# Python tools for GTMK
# Simon King, 2004
# beta version, not for redistribution
#
# ---------------------------------------------------------------


# modules
import os
import sys
import threading
from optparse import OptionParser
import gmtkPT

if __name__ == '__main__':

    # GUI-based tool for creating initial parameter files

    # future possibilities: use a .str file to actually work out what
    # parameters need to exist!

    parser = OptionParser()
    
    parser.add_option("-o", "--outputTrainableParameters", dest="param_file",
                      help="output GMTK parameter file")

    (options, args) = parser.parse_args()

    if options.param_file:
        if os.path.isfile(options.param_file):
            print "Warning: will overwrite ",options.param_file
    else:
        raise ValueError, "Must specify an output file"


    print "Starting GUI"
    gui_thread_lock=threading.Lock()        
    gui_thread=gmtkPT.Displayer(lock=gui_thread_lock,
                                tablist=['actions','makeInitialParameters'],
                                parameters={'initial_parameters_output_file' : options.param_file})


    print "Waiting for GUI to start up"
    gui_thread_lock.acquire()
    print "GUI has started up"
    gui_thread.join()
    print "Finished!"

    sys.exit(0)



