# ---------------------------------------------------------------
#
# Python tools for GTMK
# Simon King, 2004
# beta version, not for redistribution
#
# ---------------------------------------------------------------

# GMTK-format file handling

import re

def save_mlf_as_label_decision_trees(mlf,fname,mapping=None,merge=False):

    if mapping and not isinstance(mapping,dict):
        print "mapping must be a dict"
        return


    re_starts_star=re.compile(r'^(?:\*\/)?(.*)\.lab$')

    try:
        fd=open(fname,mode='w')
    except:
        print "Error opening file",fname
        #fd.close()
        return
    
    fd.write("% created by gtmkPT\n")
    fd.write("% this are decision trees giving the labels\n\n")

    if mapping:
        fd.write("% label definitions start - only used if cpp is run\n")
        fd.write(mapping['defines'])
        fd.write("% label definitions end\n\n")
                    
    fd.write(len(mlf).__repr__()+" % Number of DTs\n\n")
    c=0
    for fname,labels in mlf:
        fd.write(c.__repr__()+"\n")

        fd.write(re_starts_star.findall(fname)[0]+"_"+c.__repr__()+"\n")

        fd.write("1\n")

        # if there is a mapping, might map each label to more than one GMTK int label
        new_labels=[]
        if mapping:
            prev_m=None
            for label in labels:
                for m in mapping[label]:
                    if (not merge) or (merge and (m != prev_m)):
                        new_labels.append(m)
                    prev_m=m
        else:
            new_labels=labels


        if len(new_labels) > 1:
            fd.write("0 "+len(new_labels).__repr__())
            for r in range(len(new_labels)-1):
                fd.write(" " +r.__repr__())
                     
        else:
            fd.write("0 1")

        fd.write(" default\n")

        for label in new_labels:
            fd.write("\t-1 "+label+"\n")

        fd.write("\n")
        c=c+1
        
    fd.close()
    




def save_mlf_as_endOfUtterance_decision_trees(mlf,fname,mapping=None,merge=False):

    re_starts_star=re.compile(r'^(?:\*\/)?(.*)\.lab$')

    try:
        fd=open(fname,mode='w')
    except:
        print "Error opening file",fname
        #fd.close()
        return

    fd.write("% created by gtmkPT\n")
    fd.write("% these are decision trees for the endOfUtterance\n")
    fd.write("% parents must be: label_counter, label_transition \n\n")

    fd.write(len(mlf).__repr__()+" % Number of DTs\n\n")
    c=0
    for fname,labels in mlf:

        # if there is a mapping, might map each label to more than one GMTK int label
        new_labels=[]
        if mapping:
            prev_m=None
            for label in labels:
                for m in mapping[label]:
                    if (not merge) or (merge and (m != prev_m)):
                        new_labels.append(m)
                    prev_m=m
        else:
            new_labels=labels

        fd.write(c.__repr__()+"\n")

        fd.write('eou'+re_starts_star.findall(fname)[0]+"_"+c.__repr__()+"\n")

        fd.write("2\n")

        if len(new_labels) > 1:
            fd.write('1 2 1 default\n') # query label_transition
            fd.write('        0 2 '+(len(new_labels)-1).__repr__()+' default\n') # query label_counter
            fd.write('	          -1 1\n')
            fd.write('	          -1 0\n')
            fd.write('        -1 0\n\n')

        elif len(new_labels) == 1: # the only label IS the last label in the utterance
            fd.write('1 2 1 default\n') # just query label_transition
            fd.write('        -1 1\n')
            fd.write('        -1 0\n\n')
            
        c=c+1
        
    fd.close()
    


def save_mlf_as_label_counter_decision_trees(mlf,fname,skip_labels=[],mapping=None,merge=False):

    re_starts_star=re.compile(r'^(?:\*\/)?(.*)\.lab$')

    try:
        fd=open(fname,mode='w')
    except:
        print "Error opening file",fname
        #fd.close()
        return

    fd.write("% created by gtmkPT\n")
    fd.write("% these are decision trees for incrementing the label (e.g. word) counter\n")
    if len(skip_labels)>0:
        fd.write("% these labels can be optionally skipped; THREE parents must be: label_counter, label_transition, skip_silence\n")
        for sl in skip_labels:
            fd.write("%"+sl+"\n")
    else:
        fd.write("% no labels can be skipped; TWO parents must be: label_counter, label_transition\n")
        
    fd.write(len(mlf).__repr__()+" % Number of DTs\n\n")
    c=0
    for fname,labels in mlf:
        fd.write(c.__repr__()+"\n")
        fd.write("label_counter_"+re_starts_star.findall(fname)[0]+"_"+c.__repr__()+"\n")

        # if there is a mapping, might map each label to more than one GMTK int label
        new_labels=[]
        if mapping:
            prev_m=None
            for label in labels:
                for m in mapping[label]:
                    if (not merge) or (merge and (m != prev_m)):
                        new_labels.append(m)
                    prev_m=m
        else:
            new_labels=labels

        # assume parents are in the order: label_counter, label_transition [,skip_silence]

        if len(new_labels) < 1:
            raise ValueError, "label file with no labels:"+fname

        # 3 parents if we are skipping, otherwise 2
        if len(skip_labels)>0:
            fd.write("3\n")
        else:
            fd.write("2\n")
            
        fd.write("1 2 0 default\n") # query label_transition: either 0 or 1
        fd.write("        -1 {p0}\n") # if 0, then no transition is being made


        # a transition is being made, but is it by one or two labels?
        # (this will fail for 2 skippable labels in a row.....TO DO!?)

        # find skippable labels
        skippable=[]
        if (len(new_labels)) > 0:
            non_skippable=range(len(new_labels)-1) # exclude last label: that is handled as default case
        else:
            non_skippable=[]
        lc=1
        if len(new_labels) > 1:
            # if very last label is skippable, ignore this (nowhere to skip to)
            # also, cannot skip very first label (nowhere to skip from)
            for l in new_labels[1:len(new_labels)-2]: 
                if skip_labels.count(l) > 0:
                    # we can skip if we are in a label immediately preceding a skippable label
                    skippable.append(lc-1)
                    non_skippable.remove(lc-1)
                lc=lc+1

        skippable_string=",".join([str(skippable[x]) for x in range(len(skippable))])
        non_skippable_string=",".join([str(non_skippable[x]) for x in range(len(non_skippable))])
        final_string="default" # (len(new_labels)-1).__repr__()
        #non_skippable_string="default"

        num_branches=2 # at least the final label and the default (non-skippable) case
        if len(skippable)>0:
            num_branches=num_branches+1

            
        # query label_counter: either in skippable or non_skippable sets
        # if label is a skippable one, and skip_silence is 1, then increment word counter by 2
        fd.write("        0 "+num_branches.__repr__() \
                 +" "+skippable_string \
                 +" "+non_skippable_string \
                 +" "+final_string+"\n")

        if len(skippable)>0:
            fd.write("                2 2 0 default\n") # query skip_silence: either 0 or 1
            fd.write("                        -1 {p0+1}\n") # if 0, then no skip is being made, just a regular transition
            fd.write("                        -1 {p0+2}\n") # a skip occurs

        fd.write("                -1 {p0+1}\n") # default case, non-skippable
        fd.write("                -1 {p0}\n") # final label

        c=c+1

    fd.close()
    
