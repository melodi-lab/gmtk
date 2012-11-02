# ---------------------------------------------------------------
#
# Python tools for GTMK
# Simon King, 2004
# beta version, not for redistribution
#
# ---------------------------------------------------------------

# simple-format text file handling

import re
import os

def load_simple_text(fname):
    rl=[]
    # and return as a list of strings, one for each line in the file
    
    try:
        fd=open(fname,mode='r')
    except:
        print "Error opening file",fname
        fd.close()
        return rl

    for l in fd.readlines():
        rl.append(l)
    fd.close()
    return rl

def load_mapping(fname):

    re_split=re.compile(r'([\w\#\?\_]+)\s*')
    re_hash_defines=re.compile(r'^\#')
    re_non_blank=re.compile(r'\S')
    defines=""
    
    rval={} # a dictionary
    
    try:
        fd=open(fname,mode='r')
    except:
        print "Error opening file",fname
        fd.close()
        return rval

    for l in fd.readlines():
        if re_hash_defines.match(l):
            defines=defines+l
        elif re_non_blank.match(l):
            k = re_split.findall(l)[0]
            v = re_split.findall(l)[1:]
            rval[k]=v

    rval['defines']=defines
    fd.close()
    return rval

def load_table(fname):
    # and return as a dict
    rval={}
    re_split=re.compile(r'([\w\#\?\_\-]+)\s*')
    re_non_blank=re.compile(r'\S')

    try:
        fd=open(fname,mode='r')
    except:
        print "Error opening file",fname
        return rval

    for l in fd.readlines():
        if re_non_blank.match(l):
            k = re_split.findall(l)[0]
            v = re_split.findall(l)[1]
            rval[k]=v

    fd.close()
    return rval

def save_mlf_as_framewise_labels(outdir,mlf,framespacing,framelength,endtimes,minframes):

    re_split=re.compile(r'([\w\?\#\_\d]+)\s*')
    forced_frames=int(0)

    for fname,labels in mlf:
        bname=os.path.splitext(os.path.basename(fname))[0]
        oname=os.path.join(outdir,bname)+".ascii"
        
        try:
            fd=open(oname,mode='w')
        except:
            print "Error opening file",oname

        # assume HTK format time only, for now (TO DO: fix this)
        labs=[]
        for l in labels:
            try:
                s,e,lab=re_split.findall(l)
                st=float(s) / 10000000
                et=float(e) / 10000000
                labs.append((st,et,lab))
            except:
                print "FAILED on",l

        index=0
        t=framelength/2 # centre of first frame in HTK sig proc
        fn=0
        last_lab=""
        n_this_lab=0
        
        fet=labs[len(labs)-1][1]
        #nframes=int(fet/frameduration) - 1
        #print "frames=",nframes

        try:
            nframes=int(endtimes[bname])
        except:
            print "No end time listed for",bname

        try:
            while (fn < nframes):
                while t>labs[index][1]:
                    index=index+1
                    if index>=len(labs):
                        index=len(labs)-1
                        break
                this_lab=labs[index][2]

                if last_lab=="":
                    last_lab=this_lab # must be at start of utterance
                    
                if (this_lab != last_lab) and (n_this_lab < minframes):
                    fd.write(last_lab+"\n")
                    forced_frames=forced_frames+1
                    n_this_lab=n_this_lab+1
                else:
                    fd.write(this_lab+"\n")
                    if this_lab == last_lab:
                        n_this_lab=n_this_lab+1
                    else:
                        n_this_lab=1
                        
                    last_lab=this_lab

                t=t+framespacing
                fn=fn+1

            #print bname,fn,fet

            fd.close()
        except:
            print "Failed on",bname

    print "Forced",forced_frames,"frames"
