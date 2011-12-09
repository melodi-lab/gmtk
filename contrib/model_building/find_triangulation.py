#! /usr/bin/env python


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


# temp hacks
#gmtk_path        = '/homes/sking/projects/gmtk_dev/tksrc/'
#gmtk_triangulate = ''.join([gmtk_path,'gmtkTriangulate'])
#gmtk_tfmerge     = ''.join([gmtk_path,'gmtkTFmerge'])



def queue_boundary(options,agenda,boundary_seconds):

    kll=[]
    trifiles=[]

    # this one is the old version
    jtWeight='T'
    trifile='.'.join([options.str_file,"trifile",
                      'boundary=W',
                      'jtWeight='+jtWeight,
                      "anyTimeTriangulate=",boundary_seconds.__repr__()])
                     
    gmtk_cmd=' '.join([
        options.triangulation_cmd,
        "-strFile",options.str_file,
        "-rePartition T -reTriangulate T -findBestBoundary T",
        "-outputTriangulatedFile",trifile,
        "-triangulationHeuristic completed",
        "-boundaryHeuristic W -M 1 -S 1",
        "-seed T",
        "-numBackupFiles 0",
        "-anyTimeTriangulate",boundary_seconds.__repr__(),
        "-jtWeight T",options.cpp_arguments
        ])

    # now all the new ones
    boundaryHeuristics=['W' , 'Q', 'A',
                        'FWH', 'S', 'T',
                        'F', 'P', 'N',
                        'ST', 'SF', 'SW', 'SFW',
                        'TS', 'TF', 'TW', 'TSW',
                        'FS', 'FT', 'FW', 'FTSW']

    # try more later ?
    #'MCS', 'frontier',

    forceLeftRight=['neither','L','R']

    for jtWeight in ['T', 'F']:
        for jtwUB in ['T', 'F']:
            for b in boundaryHeuristics:
                for lr in forceLeftRight:

                    heuristic=b
                
                    trifile='.'.join([options.str_file,"trifile",
                                      'boundary='+heuristic,
                                      'jtWeight='+jtWeight,
                                      'jtwUB='+jtwUB,
                                      'forceLR='+lr])
                    
                    lrflag=''
                    if lr == 'R':
                        lrflag='-forceLeftRight R'
                    elif lr == 'L':
                        lrflag='-forceLeftRight L'

                    gmtk_cmd=' '.join([
                        options.triangulation_cmd,"-strFile",options.str_file,
                        "-rePartition T -reTriangulate T -findBestBoundary T",
                        "-boundaryHeuristic",heuristic,
                        "-jtWeight",jtWeight,
                        "-jtwUB",jtwUB,
                        lrflag,
                        "-M 1 -S 1",
                        "-seed T",
                        "-printResults T",
                        "-verbosity 20",
                        "-numBackupFiles 0",
                        "-triangulationHeuristic completed",
                        options.cpp_arguments,
                        "-outputTriangulatedFile",trifile
                        ])
                    kl=agenda.add_jobs([gmtk_cmd])
                    kll.append(kl)
                    trifiles.append(trifile)

    return kll,trifiles




def generic_queue(options,agenda,intrifile,outtrifile,heuristic,jtWeight=None,jtwUB=None):

    jtarg=''
    if jtWeight:
        jtarg=jtarg+" -jtWeight "+jtWeight
    if jtwUB:
        jtarg=jtarg+" -jtwUB "+jtwUB

    gmtk_cmd=' '.join([
        options.triangulation_cmd,"-strFile",options.str_file,
        "-rePartition F -reTriangulate T -findBestBoundary F",
        "-noReTriP F -noReTriC F -noReTriE F",
        "-inputTriangulatedFile ",intrifile,
        "-outputTriangulatedFile",outtrifile,
        "-seed T",
        "-numBackupFiles 0",
        "-triangulationHeuristic",heuristic,
        "-printResults T",
        jtarg,
        options.cpp_arguments
        ])
    return agenda.add_jobs([gmtk_cmd],[])

def queue_timing(options,agenda,trifile):
    gmtk_cmd=' '.join([options.timing_cmd,
                       "-strFile",options.str_file,
                       "-triFile",trifile,
                       "-jtFile",trifile+".jtinfo",
                       options.cpp_arguments
                       ])       
    return agenda.add_jobs([gmtk_cmd],[])

def queue_run_once_triangulations(options,agenda,trifile):
    methods=['completed', 'heuristics'] #, 'anneal']
    kll=[]

    for m in methods:
        otf='.'.join([trifile,"run_once","triangulation="+m])
        kll.append(generic_queue(options,agenda,trifile,otf,m))
    return kll

def queue_run_many_triangulations(options,agenda,trifile):
    kll=[]

    # not doing H=hint, E=entropy
    # don't yet know what these are: P, N but do them anyway

    iterations=['1', '10', '100', '1000']
    methods=['MCS', 'frontier']

    heuristics=[]
    for i in iterations:
        for m in methods:
            heuristic='-'.join([i,m])
            heuristics.append(heuristic)

    methods=[ 'S', 'T', 'F', 'W', 'X', 'P', 'N',
              'ST', 'SF', 'SW', 'SFW',
              'TS', 'TF', 'TW', 'TSW',
              'FS', 'FT', 'FW', 'FTSW']
    
    edges=['none','pre-edge-lo' ,'pre-edge-all','pre-edge-random']

    randomness=['1','2','3']

    for i in iterations:
        for e in edges:
            for m in methods:
                for r in randomness:
                    if e == 'none':
                        heuristic='-'.join([i,m,r])
                    else:
                        heuristic='-'.join([i,e,m,r])
                        heuristics.append(heuristic)

    for heuristic in heuristics:
        for jtWeight in ['T', 'F']:
            for jtwUB in ['T', 'F']:
                otf='.'.join([trifile,"triangulation="+heuristic,
                              'jtWeight='+jtWeight,
                              'jtwUB='+jtwUB])
                
                kll.append(generic_queue(options,agenda,trifile,otf,heuristic,
                                         jtWeight=jtWeight,
                                         jtwUB=jtwUB))
    return kll

def queue_run_random_triangulations(options,agenda,trifile):

    iterations=['1', '10', '100', '1000']
    methods=['R']
    attempts=range(10)

    kll=[]
    for i in iterations:
        for m in methods:
            for a in attempts:
                for jtWeight in ['T', 'F']:
                    for jtwUB in ['T', 'F']:
                        heuristic='-'.join([i,m])
                        otf='.'.join([trifile,"triangulation="+heuristic,
                                      'jtWeight='+jtWeight,
                                      'jtwUB='+jtwUB,
                                      "attempt="+a.__repr__()
                                      ])

                        kll.append(generic_queue(options,agenda,trifile,otf,heuristic,
                                                 jtWeight=jtWeight,
                                                 jtwUB=jtwUB))
    return kll




if __name__ == '__main__':

    parser = OptionParser()
    
    parser.add_option("-s", "--strfile", dest="str_file",
                      help="the structure file to be triangulated")

    parser.add_option("-c", "--cppCommandOptions", dest="cpp_arguments", default=None,
                      help="cpp options to pass to all gmtk tools")

    parser.add_option("-u", "--triangulationCommand", dest="triangulation_cmd", default="",
                      help="e.g. rexport -attr Linux gmtkTriangulate")

    parser.add_option("-g", "--timingCommand", dest="timing_cmd", default="",
                      help="e.g. rexport -attr Linux gmtkTime plus all flags except -trifile")

    parser.add_option("-b", "--boundaryTime", dest="boundary_time_string", default='20s',
                      help="how long to spend looking for a boundary")

    parser.add_option("-a", "--anyTimetriangulateTime", dest="triangulate_time_string", default='60s',
                      help="how long to run each of the anyTime triangulations for")

    # to do - make this settable for each phase separately
    parser.add_option("-m", "--maxRunTime", dest="runtime_string", default='0s',
                      help="maximum overall runtime (wall clock time); terminate any remaining jobs after this time")

    parser.add_option("-n", "--numTriangulation", dest="num_triangulations", default=0, type="int",
                      help="maximum number of triangulations to pass on to timing phase (0 means all of them)")


    parser.add_option("-p", "--maxProcesses", dest="max_processes", default=1, type="int",
                      help="maximum number of jobs to run in parallel")

    (options, args) = parser.parse_args()

    if options.str_file:
        if not os.path.isfile(options.str_file):
            raise ValueError, ''.join(["The structure file '",options.str_file,"' cannot be found"])
    else:
        raise ValueError, "Must specify a structure file"

    if options.cpp_arguments:
        options.cpp_arguments="-cppCommandOptions \\\""+options.cpp_arguments+"\\\""
        print "cpp options: ",options.cpp_arguments

    if not options.triangulation_cmd:
        print "Must give a triangulation command!"
        sys.exit(1)

    if not options.timing_cmd:
        print "Must give a timing command!"
        sys.exit(1)

    options.boundary_seconds = gmtkPT.convert_to_seconds(options.boundary_time_string)
    if options.boundary_seconds < 1:
        options.boundary_seconds=1

    options.triangulate_seconds = gmtkPT.convert_to_seconds(options.triangulate_time_string)
    if options.triangulate_seconds < 1:
        options.triangulate_seconds=1

    options.max_run_seconds = gmtkPT.convert_to_seconds(options.runtime_string)
    if options.max_run_seconds < 1:
        options.max_run_seconds=0
    
    # ================= agenda for finding best boundary ===========================
    boundary_agenda_thread_lock=threading.Lock()
    boundary_agenda_condition=threading.Condition()
    boundary_agenda_condition.acquire()
    boundary_agenda_thread=gmtkPT.Agenda(mytype='boundary',
                                         name='boundary',
                                         max_jobs=options.max_processes,
                                         lock=boundary_agenda_thread_lock,
                                         condition=boundary_agenda_condition,
                                         max_run_seconds=options.max_run_seconds)

    print "Waiting for boundary Agenda to initialise"
    boundary_agenda_thread_lock.acquire()
    print "Agenda has initialised"

    # ================= agenda for finding best triangulation ======================
    triangulation_agenda_thread_lock=threading.Lock()
    triangulation_agenda_condition=threading.Condition()
    triangulation_agenda_condition.acquire()
    triangulation_agenda_thread=gmtkPT.Agenda(mytype='triangulation',
                                              name='triangulation',
                                              max_jobs=options.max_processes,
                                              lock=triangulation_agenda_thread_lock,
                                              condition=triangulation_agenda_condition,
                                              max_run_seconds=options.max_run_seconds)

    print "Waiting for triangulation Agenda to initialise"
    triangulation_agenda_thread_lock.acquire()
    print "Agenda has initialised"


    # ================= agenda for timing best candidate triangulations ============
    timing_agenda_thread_lock=threading.Lock()
    timing_agenda_condition=threading.Condition()
    timing_agenda_condition.acquire()
    timing_agenda_thread=gmtkPT.Agenda(mytype='timing',
                                       name='timing',
                                       max_jobs=options.max_processes,
                                       lock=timing_agenda_thread_lock,
                                       condition=timing_agenda_condition,
                                       max_run_seconds=options.max_run_seconds)

    print "Waiting for timing Agenda to initialise"
    timing_agenda_thread_lock.acquire()
    print "Agenda has initialised"

    # the agendas need the lock to run
    boundary_agenda_thread_lock.release()
    triangulation_agenda_thread_lock.release()
    timing_agenda_thread_lock.release()

    print "queuing boundary jobs"
    kl,trifiles=queue_boundary(options,boundary_agenda_thread,options.boundary_seconds)
    #boundary_agenda_thread.set_end_time()

    print "Starting GUI"
    gui_thread_lock=threading.Lock()        
    gui_thread=gmtkPT.Displayer(agendas={'agenda:Boundaries':boundary_agenda_thread,
                                         'agenda:Triangulations':triangulation_agenda_thread,
                                         'agenda:Timing':timing_agenda_thread
                                         },
                                lock=gui_thread_lock,
                                tablist=['actions','updateSettings','monitor'])
    
    print "Waiting for GUI to start up"
    gui_thread_lock.acquire()
    print "GUI has started up"

   
    # signal agenda to start exporting jobs
    boundary_agenda_condition.notify()
    boundary_agenda_condition.release()
    boundary_agenda_thread.set_end_time()
    print "Boundary agenda has been signalled to go!"

    # now just wait until this Agenda has finished
    boundary_agenda_thread.join()

    # look for the top N unique most promising boundaries
    best_trifiles=boundary_agenda_thread.return_top_N_unique_trifiles()
    print "######################"
    print "from boundary algorithm, got ",len(best_trifiles),"unique trifiles :",best_trifiles
    print "######################"
    
    for tf in best_trifiles:
        queue_run_once_triangulations(options,triangulation_agenda_thread,tf)
        queue_run_many_triangulations(options,triangulation_agenda_thread,tf)
        queue_run_random_triangulations(options,triangulation_agenda_thread,tf)

    # signal agenda to start exporting jobs
    triangulation_agenda_condition.notify()
    triangulation_agenda_condition.release()
    triangulation_agenda_thread.set_end_time()
    print "Triangulation agenda has been signalled to go!"
    # now just wait until this Agenda has finished
    triangulation_agenda_thread.join()

    # look for the top N most promising triangulations
    best_trifiles=triangulation_agenda_thread.return_top_N_best_triangulations(options.num_triangulations)
    #print "got best:",best_trifiles

    for tf in best_trifiles:
        queue_timing(options,timing_agenda_thread,tf)

    # signal agenda to start exporting jobs
    timing_agenda_condition.notify()
    timing_agenda_condition.release()
    timing_agenda_thread.set_end_time()
    print "Timing agenda has been signalled to go!"
    # now just wait until this Agenda has finished
    timing_agenda_thread.join()

    gui_thread.join()
    print "Finished!"



    #agenda_thread.print_partitions()
    #print agenda_thread.return_current_best()
    sys.exit(0)
















#  def queue_PCE_triangulations(options,agenda,depends_on,trifile):

#      print "DO NOT CALL"
#      exit

#      gmtk_cmd=' '.join([
#          options.triangulation_cmd,"-strFile",options.str_file,
#          "-rePartition F -reTriangulate T -findBestBoundary F",
#          "-inputTriangulatedFile",trifile,
#          "-findBestBoundary F",
#          "-seed T",
#          "-printResults T",
#          "-numBackupFiles 0",
#          "-anyTimeTriangulate",options.triangulate_seconds.__repr__(),
#          "-jtWeight T",options.cpp_arguments,
#          "-outputTriangulatedFile"
#          ])

#      kl=agenda.add_jobs([
#          ' '.join([gmtk_cmd,'.'.join([trifile,"PCE=P"]),"-noReTriP F -noReTriC T -noReTriE T"]),
#          ' '.join([gmtk_cmd,'.'.join([trifile,"PCE=C"]),"-noReTriP T -noReTriC F -noReTriE T"]),
#          ' '.join([gmtk_cmd,'.'.join([trifile,"PCE=E"]),"-noReTriP T -noReTriC T -noReTriE F"])
#          ],depends_on)

#      gmtk_cmd=' '.join([gmtk_tfmerge,"-strFile",options.str_file,
#                         "-outputTriangulatedFile",'.'.join([trifile,"PCEmerged"]),
#                         "-Ptrifile",'.'.join([trifile,"PCE=P"]),
#                         "-Ctrifile",'.'.join([trifile,"PCE=C"]),
#                         "-Etrifile",'.'.join([trifile,"PCE=E"]),
#                         options.cpp_arguments
#                         ])


#      kl=agenda.add_jobs([gmtk_cmd],kl)

#      gmtk_cmd=' '.join([options.timing_cmd,"-strFile",options.str_file,
#                         "-triFile",'.'.join([trifile,"PCEmerged"]),
#                         options.cpp_arguments
#                         ])

#      #return agenda.add_jobs([gmtk_cmd],kl)
#      return kl

