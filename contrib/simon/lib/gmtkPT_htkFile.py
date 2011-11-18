# ---------------------------------------------------------------
#
# Python tools for GTMK
# Simon King, 2004
# beta version, not for redistribution
#
# ---------------------------------------------------------------

# HTK-format file handling

import re
import sys

def load_mlf(f,strip_times=False):
    re_is_filename=re.compile(r'^\"(.*)\".*\n$')
    re_is_delimiter=re.compile(r'^\.\n$')
    re_MLF=re.compile(r'^\#\!MLF\!\#\s*')
    re_strip_newline=re.compile(r'^(.*)\n$')
    re_split=re.compile(r'([\w\?\#\_]+)\s*')

    print "loading MLF from",f

    try:
        fd=open(f,mode='r')
        l=fd.readline()
        if not re_MLF.match(l):
            print f,"is not an MLF, first line was",l
            fd.close()
            return []
        ll=fd.readlines()
        fd.close()
    except:
        print "Error whilst reading MLF file"
        fd.close()
        return []

    # parse the MLF, putting labels for each file into a tuple
    # (filename, [label list])
    # and returning a list of such tuples
    mlf=[]
    fname=""
    labels=[]
    for l in ll:

        if re_is_filename.match(l):
            fname=re_is_filename.findall(l)[0]

        elif re_is_delimiter.match(l):
            mlf.append( (fname,labels) )
            fname=""
            labels=[]

        else:
            if strip_times:
                try:
                    s,e,lab=re_split.findall(l)
                except:
                    s=""
                    e=""
                    lab=re_strip_newline.findall(l)[0]
            else:
                lab=re_strip_newline.findall(l)[0]
                
            labels.append(lab)

    mlf.sort()
    return mlf

def save_mlf(filename=None,mlf=None):

    if not filename:
        fd=sys.stdout
    else:
        try:
            fd=open(filename,mode='w')
        except:
            print "Error opening MLF output file",filename
            return

    fd.write("#!MLF!#\n")

    for fname,labels in mlf:
        fd.write("\""+fname+"\"\n") 
        for l in labels:
            fd.write(l+"\n")
        fd.write(".\n")

    fd.close()


def mlf_mapping(mlf,mapping=None,merge=False):

    re_split=re.compile(r'([\w\?\#\_]+)\s*')

    if mapping:
        new_mlf=[]
        for fname,labels in mlf:
            
            # if there is a mapping, might map each label to more than one GMTK int label
            new_labels=[]
            prev_m=None
            for label in labels:
                # label might have start and end times
                try:
                    s,e,lab=re_split.findall(label)
                except:
                    s=""
                    e=""
                    lab=label
                    
                for m in mapping[lab]:
                    if (s=="" and e==""):
                        if (not merge) or (merge and (m != prev_m)):
                            new_labels.append(m)
                    else:
                        new_labels.append(s+" "+e+" "+m)
 
                    prev_m=m
            new_mlf.append( (fname,new_labels) )


        new_mlf.sort()
        return new_mlf

    else:
        mlf.sort()
        return mlf
