

======================================================================
======================================================================
======================================================================


			 How to compile GMTK

			  Updated: Dec 2008


Note, GMTK does not yet have a simple gnu configure style of
compilation, although that would be desirable. Therefore, to get GMTK
compiled, please follow the steps below.



1) CD to the main directory (where this file lives).

2) Edit the top level Makefile. The most important thing to set
    is the 'OPTFLAGS' flag for your appropriate architecture.
    Generic optimization is just 
     'OPTFLAGS =-g -O3 -march=pentium4 -mfpmath=sse -ffast-math'
    but you might want to use something more aggressive (but which
    if you run on the wrong architecture will produce illegal instructions).

    You might also want to change the version of gcc/g++ that you 
    use to compile.
 
3) do a 'make clean' just to make sure.

4) do a 'make depend', you will get some errors about missing
   files. Just ignore those.

5) Do a make for specific platforms. I.e., do only *one* of the following:
    make linux
    make cygwin
    make solaris
    make ibm
    make osx 

6) Do a 'make depend' once again.

7) Do 'make' and it should compile (if you are running on a 2-core
   box, do 'make -j 2' under gnumake to compile faster).

   If you are running on cygwin do 'make ANSI=' to turn off the ansi
   compatibility which is needed on that platform.

   If you wish to make static binaries (i.e., ones that do not
   depend on any shared libraries), do 'make EXLDFLAGS=-static'

   You will see an error message about 'lex.yy.c', it is safe to
   ignore it.

   You might see an error message about 'ranlib', again safe to
   ignore.

   If you wish to compile , then you can run:
    
     make MODULES=featureFileIO all 

8) All of the binaries should now live in the subdirectories:

    ./tksrc -- the main GMTK binaries, all start with 'gmtk*'

    ./featureFileIO -- a set of utilities to manipulate GMTK
                   observation files (which includes pfile, raw
                   binary, ascii, and HTK files). These can be useful
		   for preparing GMTK input observation files.
 
    ./scripts -- There are a few scripts
                 here that go through the GMTK triangulation programs
                 in a variety of ways, and they sometimes find 
                 triangulations that are significantly faster
                 than the default triangulation options.

9) Don't forget to read the documentation. It can
   be found at:

       http://ssli.ee.washington.edu/~bilmes/gmtk

10) enjoy!!

	-- Jeff Bilmes

        -- gmtk-users@ssli.ee.washington.edu
        -- gmtk-users-request@ssli.ee.washington.edu
   
======================================================================
======================================================================
