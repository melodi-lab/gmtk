# ---------------------------------------------------------------
#
# Python tools for GTMK
# Simon King, 2004
# beta version, not for redistribution
#
# ---------------------------------------------------------------

import re


def remove_empty_strings(l):
    while l.count("") > 0:
        l.remove("")
    return l


def convert_to_seconds(s):

    if not s:
        return 0
    
    s=s.lower()
    r1=re.compile(r'(\d+\s*\D+)')
    r2=re.compile(r'(\d+)\s*(\w+)')
    rweeks=re.compile(r'^w')
    rdays=re.compile(r'^d')
    rhours=re.compile(r'^h')
    rminutes=re.compile(r'^m')
    rseconds=re.compile(r'^s')

    st=r1.findall(s)
    total=0
    for i in st:
        st2=remove_empty_strings(r2.split(i))
        numbers=st2[0]; letters=st2[1]
        if rweeks.match(letters):
            total += int(numbers)*60*60*24*7
        if rdays.match(letters):
            total += int(numbers)*60*60*24
        if rhours.match(letters):
            total += int(numbers)*60*60
        if rminutes.match(letters):
            total += int(numbers)*60
        if rseconds.match(letters):
            total += int(numbers)

    return total # seconds

# some useful regexs
class Regex:
    def __init__(self):
        self.remove_newline=re.compile(r'(.*)')
        self.remove_path=re.compile(r'.*\/(\w+)')
        self.trifile_name=re.compile(r'outputTriangulatedFile\s+([\w*\.*\/*\=*\-*]+)') # matches gmtkTriangulate args
        self.trifile_name_alt=re.compile(r'triFile\s+([\w*\.*\/*\=*\-*]+)') # matches gmtkTime args
        self.is_gmtkTime=re.compile(r'.*gmtkTime.*')
        self.is_gmtkTriangulate=re.compile(r'.*gmtkTriangulate.*')
        self.is_gmtkTriangulate_find_boundary=re.compile(r'.*gmtkTriangulate.*boundaryHeuristic.*')
        self.partitions=re.compile(r'([0-9]+) total partitions')
        self.jtWeights=re.compile(r'jt_weight = ([0-9\.]+)')

        self.which_interface=re.compile(r'Using (\w+) interface to define partitions')
        self.best_interface_line=re.compile(r'.*Best interface nodes include.*')
        self.best_interfaces=re.compile(r'(\S+\(\d+\))')
