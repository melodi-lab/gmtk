
EXCCFLAGS=
EXCXXFLAGS=
EXLDFLAGS =  

OPTFLAGS = -DHAVE_NONSTANDARD_ARITHMETIC

ANSI=-ansi
PEDANTIC=-pedantic
WALL=-Wall

# Build libmiscSupport.a, but don't install it in $(libdir)
noinst_LIBRARIES = libmiscSupport.a

# The source files libmiscSupport.a depends on
libmiscSupport_a_SOURCES = error.h error.cc bp_range.h bp_range.cc range.h range.cc arguments.h arguments.cc logp.h logp.cc rand.h rand.cc sArray.h sArray.cc cArray.h cArray.cc general.h general.cc lineqsolve.h lineqsolve.cc debug.h debug.cc fileParser.h fileParser.cc hash_abstract.h hash_abstract.cc mul_mdmd_md.h mul_mdmd_md.c mul_mdmd_md_l0g.h mul_mdmd_md_l0g.c mul_mdmd_md_l0nf.h mul_mdmd_md_l0nf.c eig.h eig.c bin_search.h bin_search.cc lms_filter.h lms_filter.cc rls_filter.h rls_filter.cc

# Compilation flags for miscSupport library objects (user's CXXFLAGS will
# be appended, defaults to "-g -O2" for gcc)
libmiscSupport_a_CXXFLAGS = $(OPTFLAGS) $(EXCXXFLAGS) $(WALL) $(ANSI) $(PEDANTIC)
libmiscSupport_a_CCFLAGS = $(OPTFLAGS) $(EXCCFLAGS) $(WALL) $(ANSI) $(PEDANTIC)

# These are test programs, so they only need to be built for make check.
# They are not installed in $(bindir)

# testError testDebug testBPRange 
check_PROGRAMS = testArguments testLogp
# testRand testsArray testcArray fileParserTest testHashDrivers testHashMapList testHashTree testVHashSet testVHashMap testSHashSet testSHashMap testVSHashMap

# Complilation flags for test programs (user's CXXFLAGS will be appended)
AM_CXXFLAGS = $(WALL) -DMAIN
# Linker flags for test programs (user's LDFLAGS will be appended)
AM_LDFLAGS = $(EXLDFLAGS)

#testError_SOURCES = error.cc
#testError_LDADD = $(objs)

#testDebug_SOURCES = debug.cc debug.h
#testDebug_LDADD = $(objs)

#testBPRange_SOURCES = bp_range.cc 
#testBPRange_LDADD = error.o $(objs)

testArguments_SOURCES = arguments.cc
testArguments_LDADD = libmiscSupport.a

testLogp_SOURCES = logp.cc
testLogp_LDADD = libmiscSupport.a

testRand: rand.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ rand.cc error.o $(objs)

testsArray: sArray.cc sArray.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ sArray.cc error.o $(objs)

testcArray: cArray.cc
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ cArray.cc error.o $(objs)


fileParserTest: fileParser.cc fileParser.h
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ fileParser.cc $(objs)

testHashDrivers: testHashMapList testHashTree testVHashSet testVHashMap testSHashSet testSHashMap testVSHashMap

testHashMapList: hash_map_list.cc hash_abstract.o hash_map_list.h hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ hash_map_list.cc error.o hash_abstract.o general.o

testHashTree: hash_tree.h hash_tree.cc hash_map_list.cc hash_abstract.o hash_map_list.h hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ hash_tree.cc error.o hash_abstract.o general.o

testVHashSet: vhash_set.h vhash_set.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ vhash_set.cc error.o hash_abstract.o general.o

testVHashMap: vhash_map.h vhash_map.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN  -o $@ vhash_map.cc error.o hash_abstract.o general.o

testSHashSet: shash_set.h shash_set.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN $(HASH) -o $@ shash_set.cc error.o hash_abstract.cc general.o 

testSHashMap: shash_map.h shash_map.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN $(HASH) -o $@ shash_map.cc error.o hash_abstract.cc general.o

testVSHashMap: vshash_map.h vshash_map.cc hash_abstract.o hash_abstract.cc hash_abstract.h general.o $(objs)
	$(compile.cc.link) -Wall -g -DMAIN $(HASH) -o $@ vshash_map.cc error.o hash_abstract.cc general.o

testRLS: rls_filter.cc rls_filter.h  adaptive_filter.h libmiscSupport.a
	$(compile.cc.link) -Wall -g -DMAIN -o $@ rls_filter.cc libmiscSupport.a

testLMS: lms_filter.cc lms_filter.h  adaptive_filter.h libmiscSupport.a
	$(compile.cc.link) -Wall -g -DMAIN -o $@ lms_filter.cc libmiscSupport.a

