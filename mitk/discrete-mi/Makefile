#CC	=	/g/ssli/commercial/pure/purify-5.1-solaris2/purify gcc
#C++	=	/g/ssli/commercial/pure/purify-5.1-solaris2/purify g++

CC	=	gcc
C++	=	g++

ROOT=../..

INCLUDE	=	-I. -I$(ROOT)/miscSupport -I$(ROOT)/featureFileIO -I$(ROOT)/IEEEFloatingpoint
LIBDIR	=	-L. -L$(ROOT)/miscSupport -L$(ROOT)/featureFileIO
LDLIBS	= 	-lgmtkio -lmiscSupport -lm 

DEBUG	=	-g
OPTIMIZE =	-O3

CFLAGS   = 	$(DEBUG) $(OPTIMIZE) -Wall
CPPFLAGS =	$(DEBUG) $(OPTIMIZE) -Wall

######################################################################
.SUFFIXES: .c .cc .cpp

COMPILE.c=$(CC) $(INCLUDE) $(CFLAGS) -c
COMPILE.cpp=$(C++) $(CPPFLAGS) $(INCLUDE) -c
LINK = $(C++) $(CPPFLAGS) $(LIBDIR)
######################################################################

PROJECT1 = discrete-mi

all: $(PROJECT1)

OBJS1 = mixNormalCollection.o discrete_mi.o mixNormal.o readRange.o rangeSetCollection.o

$(PROJECT1): $(OBJS1)
	$(LINK) -o $@ $(OBJS1) $(LDLIBS)


######################################################################
.c:
	$(COMPILE.c) -c $<

.c.o:
	$(COMPILE.c) -c $<

.cpp.o:
	$(COMPILE.cpp) -c $<

.cc.o:
	$(COMPILE.cpp) -c $<
######################################################################

check-in:
	ci -u *.h *.H *.cc Makefile

check-out:
	co -l *.h *.H *.cc Makefile

######################################################################

clean:
	rm -f *.o core









