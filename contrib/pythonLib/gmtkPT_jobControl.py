# ---------------------------------------------------------------
#
# Python tools for GTMK
# Simon King, 2004
# beta version, not for redistribution
#
# ---------------------------------------------------------------

# modules
import os
import time
import re
import threading
import popen2
import signal
import pwd
import Tkinter

import gmtkPT

class Status:
    "pseudo enumerated type"
    notStarted,sentToPmake,running,completed,failed = range(5)
    text=["not started","sent to pmake","running","complete","failed"]
    text_brief=["N","P","R","C","F"]
    
class Job:
    "used by agenda"
    command="None"
    status=Status.notStarted
    pid=0
    returnValue=0
    output=""
    partitions=0
    sum_jtWeight=None
    node="Unknown"
    depends_on=None # which other jobs must finish before starting this one (list of keys into Agenda.jobs)

    def __init__(self,c,depends_on=[]):

        self.command=c
        self.depends_on=depends_on
        
    def printme(self,key):
        print "JOB",key,":",self.command
        print " status:",Status.text[self.status]
        print " node:  ",self.node
        print " pid:   ",self.pid
        print " rval:  ",self.returnValue
        print " output:",self.output
        print " parts: ",self.partitions
        print "depends:",self.depends_on

    def printme_brief(self,key):
        # to do : rewrite to call return_brief_summary
        oput=""
        if len(self.output) > 0:
            oput=gmtkPT.regex.remove_newline.findall(self.output[len(self.output)-1])[0]

        jstr=""
        if len(self.depends_on) > 0:
            jstr=self.depends_on[0].__repr__()
            for i in self.depends_on[1:]:
                jstr=jstr+","+i.__repr__()

        pstr=""
        if self.partitions > 0:
            pstr = "["+self.partitions.__repr__()+"]"
        try:
            cmd=gmtkPT.regex.remove_path.findall(self.command.split()[0])[0]
        except IndexError:
            cmd="No command yet"
        return key.__repr__() \
               +" "+Status.text_brief[self.status] \
               +" ("+jstr+")" \
               +pstr \
               +" "+self.node \
               +" "+cmd \
               +" "+oput


    def return_brief_summary(self,key):

        oput=""
        if len(self.output) > 0:
            oput=gmtkPT.regex.remove_newline.findall(self.output[len(self.output)-1])[0]

        jstr=""
        if len(self.depends_on) > 0:
            jstr=self.depends_on[0].__repr__()
            for i in self.depends_on[1:]:
                jstr=jstr+","+i.__repr__()

        pstr=""
        if self.partitions > 0:
            pstr = "["+self.partitions.__repr__()+"]"
        elif self.sum_jtWeight != None:
            pstr = "("+self.sum_jtWeight.__repr__()+")"
        else:
            pstr = "-"

        return (key.__repr__(), \
                Status.text_brief[self.status], \
                jstr, \
                pstr, \
                self.node, \
                self.command, \
                oput)

    #gmtkPT.regex.remove_path.findall(self.command.split()[0])[0], \

    def save_stats(self,append=False):
        # find the trifile name
        tfile=gmtkPT.regex.trifile_name.findall(self.command)
        if len(tfile) == 0:
            tfile=gmtkPT.regex.trifile_name_alt.findall(self.command)
        if len(tfile) == 0:
            print "Cannot save - no trifile to be found from"
            self.printme('unknown')
            exit
        else:
            trifilename=tfile[0]
            statsfilename=trifilename+".stats"
            #print "saving to",statsfilename

            if append:
                fd=open(statsfilename,'a')
                fd.write("(Appending to file)\n")
            else:
                fd=open(statsfilename,'w')
                fd.write("Stats for trifile ")
                fd.write(trifilename)

            fd.write("\n\nCreated by ")
            fd.write(self.command)
            fd.write("\n\n")

            fd.write("returnValue = ")
            fd.write(self.returnValue.__repr__())
            fd.write("\n")

            # TO DO: don't sleep if we are in the process of exiting 
            time.sleep(1)
            fd.writelines(self.output)
            fd.write("\n")
            time.sleep(1)

            fd.write("status: ")
            fd.write(Status.text[self.status])
            fd.write("\n")

            if self.partitions > 0:
                fd.write("partitions: ")
                fd.write(self.partitions.__repr__())
                fd.write("\n")

            if self.sum_jtWeight:
                fd.write("sum_jtWeight (P+C+E): ")
                fd.write(self.sum_jtWeight.__repr__())
                fd.write("\n")

            fd.close()
            #time.sleep(1)



class myThread(threading.Thread):
    re1=None
    re2=None
    fo=None
    pobj=None
    job=None
    def __init__(self, target=None, args=(), name=None):
        self.re1=re.compile(r'.*exported to\s*(\S*)\.ee\.washington\.edu \s*.*')
        self.re2=re.compile(r'.*\w+.*')
        target=self.run_rexport
        threading.Thread.__init__(self, target=target, name=name, args=args)
        self.setDaemon(False)
        self.start()

    def __del__(self):
        print "This is del boy"

    def die(self):
        if self.pobj != None:
            print "kill -"+signal.SIGKILL.__repr__(),self.pobj.pid
            try:
                os.kill(self.pobj.pid,signal.SIGKILL)
            except:
                print "didn't manage to kill process",self.pobj.pid
        else:
            print "no PID here"

    def print_info(self):
        # print which node we exported to
        if self.job:
            if len(self.job.output) > 0:
                return self.job.node+" "+self.job.output[len(self.job.output) - 1]
            else:
                return self.job.node+" no output yet"
        else:
            return "no job"
        
    def run_rexport(self,c,j,a):
        self.job=j
        #print "set job to",self.job,j
        j.status=Status.running
        j.output=[]
        self.pobj=popen2.Popen4(c)
        
        while self.pobj.poll() < 0:
            
            line="dummy"
            while line != "":
                try:
                    line = self.pobj.fromchild.readline()
                    m=self.re1.findall(line)
                    if len(m) > 0:
                        self.job.node=m[0]
                        #print "thread read",line
                    if self.re2.match(line):
                        j.output.append(line)
                    #else:
                    #    print "no match line:",line
                except:
                    line = ""
                    
            # sleep (this thread only)
            time.sleep(5) # maybe make this configurable ?

        if self.pobj.poll() != 9:
            j.status=Status.completed
            lines=self.pobj.fromchild.readlines() # this will hang for killed processes??
            for l in lines:
                time.sleep(0.05)
                j.output.append(l)
        else:
            j.status=Status.failed
            j.output.append("No more output: process was killed")

        self.pobj.tochild.close()
        self.pobj.fromchild.close()
        j.returnValue=self.pobj.poll()
        # now grab any interesting info from the output
        if gmtkPT.regex.is_gmtkTime.match(c):
            m=gmtkPT.regex.partitions.findall(''.join(j.output))
            if len(m) > 0:
                j.partitions=int(m[0])

        elif gmtkPT.regex.is_gmtkTriangulate_find_boundary.match(c):
            j.sum_jtWeight=-1 # not keeping jt weight of existing trianguklation
            
        elif gmtkPT.regex.is_gmtkTriangulate.match(c):
            m=gmtkPT.regex.jtWeights.findall(' '.join(j.output))
            if len(m) > 0:
                j.jtWeights=' , '.join(m)
                j.sum_jtWeight=float(0.0)
                for mi in m:
                    j.sum_jtWeight=j.sum_jtWeight+float(mi)

        # delete popen4 object
        self.pobj=None

        # and save this info in a file
        self.job.save_stats(append=False)

class Agenda(threading.Thread):
    "a list of jobs to be done, in progress, and completed"
    next_key=0
    jobs=dict()
    sleep_time=1 # to do: make this user settable
    end_time=None
    start_time=None
    terminate_now=False
    
    def __init__(self,
                 mytype='',
                 max_jobs=1,
                 lock=None, condition=None,
                 max_run_seconds=0,
                 target=None, name='gmtkPT_Agenda'):

        #self.rexport=''.join(["rexport ",rexport])
        #self.rexport_timing=''.join(["rexport ",rexport_timing])

        #print "rexport is",rexport

        #print "rexport_timing is",rexport_timing

        self.mytype=mytype
        self.max_run_seconds=max_run_seconds

        print "Created agenda type:",self.mytype

        self.thread_lock=lock
        if self.thread_lock:
            self.thread_lock.acquire()
        else:
            print "Warning: Agenda starting without thread locking"

        self.condition=condition
        if not self.condition:
            print "Warning: Agenda starting without condition signalling"

        if not target:
            target=self.run

        if max_jobs < 1:
            print "Can't have fewer than 1 job, setting max_jobs to 1"
            max_jobs=1

        self.max_threads=max_jobs
        self.jobs={}

        # start this thread
        threading.Thread.__init__(self, target=target, name=name, args=())
        self.setDaemon(False)
        self.start()

        # signal to calling thread that we are initialised and ready to go
        if self.thread_lock:
            self.thread_lock.release()


    def set_end_time(self):
        # work out the wall-clock end time
        if self.max_run_seconds > 0:
            self.end_time=time.time() + self.max_run_seconds
            print "Local time now is ",time.asctime(time.localtime(time.time()))
            print "Agenda",self.getName(),"will terminate no later than",time.asctime(time.localtime(self.end_time))
        else:
            self.end_time=None
            print "Will terminate only when all jobs are completed"

        self.start_time=time.time()

    def print_info(self):
        print "Agenda thread:"
        print " name",self.getName()
        print " jobs",len(self.jobs)
        print " end time",self.end_time
        
    def add_job(self,c,depends_on=[]):
        #print "agenda",self.getName(),"adding job",c
        self.jobs[self.next_key]=Job(c,depends_on)
        k=self.next_key
        self.next_key+=1
        return k

    def add_jobs(self,cl,depends_on=[]):
        kl=[]
        for c in cl:
            k=self.add_job(c,depends_on)
            kl.append(k)
        return kl

    def printme(self):
        print "---AGENDA---"
        for key in self.jobs:
            self.jobs[key].printme(key)

    def printme_brief(self):
        print "---AGENDA---"
        for key in self.jobs:
            print self.jobs[key].printme_brief(key)

    def return_brief_summary(self,order=None,list=None,
                             num_total=0,num_completed=0):
        l=[]

        if order and list:
            pos_in_list=[None for x in range(len(order))] # list of "None"s of correct length
            index=0
            for i in order:
                pos_in_list[int(i)]=index
                index=index+1

            for key in self.jobs:
                num_total=num_total+1
                if self.jobs[key].status == Status.completed:
                    num_completed=num_completed+1
                s=self.jobs[key].return_brief_summary(key)
                l.append(s)
                try:
                    pos=pos_in_list[key] # throws exception for key out of range
                    if pos != None:
                        list.delete(pos,pos)
                        list.insert(pos,s)
                    else:
                        list.insert(Tkinter.END,s)

                except IndexError:
                    print "New agenda item:",key,s
                    list.insert(Tkinter.END,s)

            return (num_completed,num_total)


        elif list:
            print "This is the one"
            list.delete(0,Tkinter.END)
            for key in self.jobs:
                num_total=num_total+1
                if self.jobs[key].status == Status.completed:
                    num_completed=num_completed+1
                s=self.jobs[key].return_brief_summary(key)
                l.append(s)
                list.insert(Tkinter.END,s)
            return (num_completed,num_total)
            
        else:
            l=[]
            for key in self.jobs:
                num_total=num_total+1
                if self.jobs[key].status == Status.completed:
                    num_completed=num_completed+1
                l.append(self.jobs[key].return_brief_summary(key))
            return (num_completed,num_total)

                
    def print_partitions(self):
        print "---Partitions---"
        for key in self.jobs:
            if self.jobs[key].partitions > 0:
                print key.__repr__(),self.jobs[key].return_brief_summary(key)

    def remaining_jobs(self):
        c=0
        for key in self.jobs:
            if self.jobs[key].status < Status.completed:
                c=c+1
        return c

    def kill_all_processes(self):
        for t in threading.enumerate():
            if t.getName() != "MainThread" :
                print " killing thread ",t.getName()
                t.die()

    def export_jobs_via_threads(self,job_key_list):
        re1=re.compile(r'gmtkPT_job.*')
        re2=re.compile(r'.*gmtkTime.*') # could replace with standard regex

        active_threads=0
        for t in threading.enumerate():
            if re1.match(t.getName()):
                active_threads=active_threads+1

        for job_key in job_key_list:
            if active_threads < self.max_threads:
                if self.jobs[job_key].status == Status.notStarted:
                    
                    depends_on_complete=True
                    for k in self.jobs[job_key].depends_on:
                        if self.jobs[k].status != Status.completed:
                            depends_on_complete=False
                            
                    if depends_on_complete:
                        print "starting thread for job",job_key
                        
                        c=' '.join([self.jobs[job_key].command,'2>&1'])

                        print c
                        t=myThread(name='gmtkPT_job'+job_key.__repr__(),args=[c,self.jobs[job_key],self])
                        active_threads=active_threads+1

                        time.sleep(self.sleep_time) 


    def watch_threads(self):
        print "active threads:",threading.activeCount()
        for t in threading.enumerate():
            if t.getName() != "MainThread":
                print " ",t.getName(),t.print_info()
                #self.find_current_best()

    def return_current_best(self):
        r=""
        # look for best number of partitions or best sum_jtWeight
        bp=0
        for key in self.jobs:
            if self.jobs[key].partitions > bp:
                bp=self.jobs[key].partitions
                bk=key
        bjt=9999
        for key in self.jobs:
            if self.jobs[key].sum_jtWeight and self.jobs[key].sum_jtWeight < bjt:
                bjt=self.jobs[key].sum_jtWeight
                bjtk=key

        if bp>0:
            ckeys=self.find_calling_job(bk)
            r=r+"Best partitions="+bp.__repr__()+"\n"
            r=r+" from job "+bk.__repr__()+"\n"
            if len(ckeys) > 0:
                r=r+" called by job(s)\n"
                for k in ckeys:
                    r=r+" "+k.__repr__()+" "+self.jobs[k].command+"\n"

        elif bjt<9999 and bjt>0:
            ckeys=self.find_calling_job(bjtk)
            r=r+"Best sum_jtWeight="+bjt.__repr__()+"\n"
            r=r+" from job "+bjtk.__repr__()+"\n"
            if len(ckeys) > 0:
                r=r+" called by job(s)\n"
                for k in ckeys:
                    r=r+" "+k.__repr__()+" "+self.jobs[k].command+"\n"


        else:
            r="No current best....yet"

        return r

    def unique(self,interfaces):
        # interfaces must be a dict of "interface" keyed by "key"
        unique_interfaces={} # a dict of "key" keyed by "interface"
        for key in interfaces:
            if not unique_interfaces.has_key(interfaces[key]):
                unique_interfaces[interfaces[key]]=key
        return unique_interfaces

    def dict_to_list(self,d,reverse=False):
        # convert dict to list of 2-tuples (key,value)
        l=[]
        for k in d:
            if reverse:
                l.append((d,k))
            else:
                l.append((k,d))
        return l

    def grab_interface(self, output):
        # "output" is a list of strings, each being one line of gmtkTriangulate output
        # return value is a list of RVs making up the interface
        
        # output will look something like:
        # Finding BEST LEFT interface
        # New Best Triangulation Found:          1-completed 2.528917  
        # etc
        # Best interface nodes include: word(2) wordPosition(2)
        # Finding BEST RIGHT interface
        # New Best Triangulation Found:          1-completed 2.528917  
        # etc
        # Best interface nodes include: word(2) wordPosition(2) phoneTransition(2) wordTransition(2)
        # etc
        # Using left interface to define partitions

        # so, first look for "Using [left|right] interface to define partitions"
        oneline=''.join(output)

        m=gmtkPT.regex.which_interface.findall(oneline)
        if len(m) != 1:
            #print "Cannot work out which interface was chosen:",m
            return None
        else:
            which_interface=m[0]
            if which_interface == "left":
                which_interface = "LEFT"
            elif which_interface == "right":
                which_interface = "RIGHT"
            else:
                print "Unexpected value for which interface:",which_interface
                return None

        # now, look for "Finding BEST which_interface interface"
        re1=re.compile(r'Finding BEST '+which_interface+' interface')
        flag=False
        i=None
        for l in output:
            if re1.match(l):
                flag=True
            if flag and gmtkPT.regex.best_interface_line.match(l):
                i=gmtkPT.regex.best_interfaces.findall(l)
                break

        return i

    def list_to_string(self,l):
        s=''
        for i in l:
            s=s+i+' '
        return s
    
    def return_top_N_unique_trifiles(self,N=0):
        # for now, return ALL unique files (no proper way of ranking them)
        if self.mytype != 'boundary':
            print "Can't return_top_N_unique_trifiles from this agenda",self.getName()
            return []

        interfaces={}
        for key in self.jobs:
            interface=self.grab_interface(self.jobs[key].output)
            if interface:
                interface.sort()
                interfaces[key]=self.list_to_string(interface)

        #print "all interfaces:"
        #print interfaces
        #print
        
        interfaces=self.unique(interfaces)
        # convert to list, if need to sort
        #sort_by_size(interfaces) # so we process smaller interfaces first - might be faster??

        #print "unique interfaces:"
        #print interfaces
        #print
        

        #if N==0 or N>len(interfaces):
        #    N=len(interfaces)

        tf=[]
        re1=re.compile(r'\-outputTriangulatedFile (\S+)')
        for key in interfaces:
            # append the trifile corresponding to this key
            t=re1.findall(self.jobs[interfaces[key]].command)[0]
            #print "this trifile will be passed to next phase:",t
            tf.append(t)
        return tf


    def return_top_N_best_triangulations(self,N=0):
        if self.mytype != 'triangulation':
            print "Can't return_top_N_best_triangulations from this agenda",self.getName()
            return []


        jtw_list=[]
        for key in self.jobs:
            if self.jobs[key].sum_jtWeight:
                jtw_list.append((self.jobs[key].sum_jtWeight,key))

        print "Got a list of",len(jtw_list),"jt weights"
        jtw_list.sort()
        print jtw_list

        if N==0 or N>len(jtw_list):
            N=len(jtw_list)


        if N==0:
            #no jt weights at all
            return []
            
        re1=re.compile(r'\-outputTriangulatedFile (\S+)')
        tf=[]
        threshold=jtw_list[N-1][0]

        print "threshold",threshold

        for i in range(N): 
            key=jtw_list[i][1]
            print "appending",i,key,jtw_list[i][0]
            tf.append(re1.findall(self.jobs[key].command)[0])

        print "Top N:",

        print "returning",N,"jobs of a possible",len(jtw_list)
        print "best/worst jtWeights are:",jtw_list[0][0],threshold
        return tf



    def find_calling_job(self,key):
        if len(self.jobs[key].depends_on) == 0:
            print "Cannot determine calling job for:",key
            #self.jobs[key].printme(key)
            return []
        else:
            return self.jobs[key].depends_on


    def run(self):
        re1=re.compile(r'gmtkPT_job.*')

        # wait for the signal
        print "Agenda is waiting to go"
        self.condition.acquire()
        print "Agenda is running"

        # the main loop for this class: repeatedly exports jobs until all are finished
        active_threads=0
        while (self.remaining_jobs() > 0) or (active_threads > 0):
            #print "exporting...."

            # need the lock whilst exporting so we don't export while the GUI is exiting
            if self.thread_lock:
                self.thread_lock.acquire()

            if self.terminate_now:
                print "Agenda",self.getName(),"terminated"
                if self.thread_lock:
                    self.thread_lock.release()
                return

            if self.end_time and (time.time() > self.end_time):
                print "Agenda",self.getName(),"ran out of time at",time.asctime(time.localtime(time.time()))

                # kill the threads with jobs in my list
                for job_key in self.jobs.keys():
                    name='gmtkPT_job'+job_key.__repr__()
                    for t in threading.enumerate():
                        if t.getName() == name :
                            print " killing remaining thread ",t.getName()
                            t.job.output.append("Aborted: ran out of time")
                            t.die()


                self.terminate_now=True
                if self.thread_lock:
                    self.thread_lock.release()
                return

            time.sleep(self.sleep_time)

            self.export_jobs_via_threads(self.jobs.keys())

            if self.thread_lock:
                self.thread_lock.release()

            time.sleep(self.sleep_time)
            active_threads=0
            for t in threading.enumerate():
                time.sleep(self.sleep_time / 10) # mini-sleep
                if re1.match(t.getName()):
                    active_threads=active_threads+1
                    #print "there are now",active_threads,"job threads"
        print "Agenda completed"
