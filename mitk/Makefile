CC	=	gcc
C++	=	g++
#C++	=	/uns/bin/g++

INCLUDES= -I.  -I../miscSupport  -I../featureFileIO -I../IEEEFloatingpoint
LIBDIR= -L. -L../miscSupport -L../featureFileIO -L../IEEEFloatingpoint
LIBS=-lgmtkio -lmiscSupport -lIEEEsupport -lm

OPTFLAGS =	-O3

CFLAGS	=	-g $(OPTFLAGS) $(EXCFLAGS) -Wall -Wno-deprecated -ansi -pedantic
CPPFLAGS =	-g $(OPTFLAGS) $(EXCFLAGS) -Wall -Wno-deprecated -ansi -pedantic

ifdef DEBUG
CPPFLAGS := $(CPPFLAGS) -DDEBUG='$(DEBUG)'
CFLAGS := $(CFLAGS) -DDEBUG='$(DEBUG)'
endif

######################################################################
.SUFFIXES: .c .cc .cpp

COMPILE.c=$(CC) $(INCLUDES) $(CFLAGS)
COMPILE.cpp=$(C++) $(CPPFLAGS) $(INCLUDES)
LINK = $(C++) $(CPPFLAGS) $(LIBDIR)

######################################################################

MULVMI  = multivariate-mi
BIVMI   = bivariate-mi

all: $(MULVMI) $(BIVMI) 

# Add the below if you need to test them
#test-matrix-ops mvnrnd testReadRange

mvnrnd: mvnrnd.o matrix-ops.o
	$(LINK) -o $@ mvnrnd.o matrix-ops.o $(LIBS)


testReadRange: testReadRange.o readRange.o rangeSetCollection.o
	$(LINK) -o $@ testReadRange.o readRange.o rangeSetCollection.o $(LIBS)


OBJS = mixNormalCollection.o multivariate_mi.o mixNormal.o kmeans.o readRange.o rangeSetCollection.o partition.o computeMI.o computeMI-data.o data-points.o matrix-ops.o tests.o io.o generate-synthetic-data.o compute-entropy.o

BIVMI_OBJS = bivariate_mi.o MixBiNormal.o MixBiNormal_chelp.o


SRCS =  $(OBJS:.o=.cc)

$(MULVMI): $(OBJS)
	$(LINK) -o $@ $(OBJS) $(LIBS)


$(BIVMI): $(BIVMI_OBJS)
	$(LINK) -o $@ $(BIVMI_OBJS) $(LIBS)

test-matrix-ops: test-matrix-ops.o matrix-ops.o
	$(LINK) -o $@ test-matrix-ops.o matrix-ops.o $(LIBS)

#
# You shouldn't have to change anything below here.
#
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

clean:
	rm -f *.o core $(MULVMI) $(BIVMI)

######################################################################

depend:
	${CXX} ${INCLUDES} -MM ${CXXFLAGS} ${SRCS} > depends.make

depends: depend

##################
####  DEPENDS ####
include depends.make
