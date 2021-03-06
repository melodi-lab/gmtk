

<img src="/images/gmtk_logo.png" align="center" width="75%">


The Graphical Models Toolkit (GMTK) is an open source, publicly available toolkit for rapidly prototyping statistical models using dynamic graphical models (DGMs) and dynamic Bayesian networks (DBNs). GMTK can be used for applications and research in speech and language processing, bioinformatics, activity recognition, and any time series application. GMTK has many features, including exact and approximate inference; a large variety of built-in factors including dense, sparse, and deterministic conditional probability tables, native support for ARPA backoff-based factors and factored language models, parameter sharing, gamma and beta distributions, dense and sparse Gaussian factors, heterogeneous mixtures, deep neural network factors, and time-inhomogeneous trellis factors; arbitrary order embedded Markov chains; a GUI-based graph viewer; flexible feature-file support and processing tools (supporting pfiles, HTK files, ASCII/binary, and HDF5 files); and both offline and streaming online inference methods that can be used for both parameter learning and prediction. More information is available in the documentation. All in all, GMTK offers a flexible, concise, and expressive probabilistic modeling framework with which one may rapidly specify a vast collection of temporal statistical models.

GMTK was developed by Jeff Bilmes, Richard Rogers, and a number of other individuals. Please see the PDF documentation for complete details and acknowledgments. Work on GMTK was supported by NIH award U01 HG009395 and by the CONIX Research Center, one of six centers in JUMP, a Semiconductor Research Corporation (SRC) program sponsored by DARPA. Support was also provided by the TerraSwarm Research Center, one of six centers administered by the STARnet phase of the Focus Center Research Program (FCRP) a Semiconductor Research Corporation program sponsored by MARCO and DARPA, the National Science Foundation grants CNS-0855230, IIS-0905341, IIS-0093430, IIS-0434720, and IIS-0326382, DARPA's ASSIST Program (contract number NBCH-C-05-0137), NIH awards R01 GM096306 and P41 GM103533, an ONR MURI grant (No. N000140510388), and generous gifts by Microsoft Research, IBM, the Intel corporation, and Google.


## Documentation

Documentation for GMTK, and on dynamic graphical models in general, [is available in this PDF file](./documentation.pdf). While the documentation is not complete, it is at this point over 600 pages and contains some hopefully useful information.

#  How to compile GMTK

GMTK uses a simple GNU Autotools build process. The familiar

```
./configure && make && make install
```

should work on the target platforms (Linux, Mac OS X, and Windows/Cygwin).
Other POSIX-like platforms might work too. By default, the GMTK programs
will be installed to /usr/local/bin. You can change that with the --prefix
option to the configure script. See

```
./configure --help
```

for full information on the available configure options.


The program `gmtkDMLPtrain` requires a compiler with ISO C++ 2011 support. The configure
script will enable support for GCC 4.5 or later. This may also work for
other compilers that enable ISO C++ 2011 support with `-std=c++0x` or `-std=c++11`.
For compilers for which that doesn't work, you may need to specify the 
required command line arguments with something like

```
./configure CXXFLAGS='<arguments to enable C++ 2011>' ...
```

If you build GMTK with a compiler that does not support ISO C++ 2011,
the `gmtkDMLPtrain` program will not be compiled. You will also want to 
use a good CBLAS implementation like MKL or ATLAS to get the best
performance with `gmtkDMLPtrain`. See the relevant options in `./configure --help`

Also note that on some platforms, ABI incompatibilities can cause
programs that mix ISO C++ 2011 and non-2011 C++ to fail (they may
compile & link successfully, but crash or produce incorrect output).
`gmtkDMLPtrain` may not work properly on such platforms.


`gmtkViz` (which is a graphical user interface for taking an existing graph structure specified using a .str file, and quickly producing a nice layout of that graph, and a resulting .eps/.pdf file for, say, inclusion in a paper) is only built if a suitable [wxWidgets](http://www.wxwidgets.org/) installation version 3.0 or later is detected at configure time. The "Print to EPS file" option under the file menu does not work with wxWidgets 2.9.2 or earlier due to bugs in wxWidgets. You can use the native OS print-to-file driver to produce a PDF file, but on Linux it seems to generate bitmap rather than vector output. Note that wxWidgets is the only external dependence needed to compile any of GMTK, nor is it required to use GMTK.

On most platforms, gmtkViz will be built automatically if a suitable 
wx-config is available under your PATH environment variable. You can
also specify how to find wxWidgets using the `--with-wx` configure options.
See `configure --help` for more details. If a suitable wxWidgets installation
is not found, `gmtkViz` will be skipped but the rest of GMTK will build and
work fine.

To build `gmtkViz` on Mac OSX we recommend installing wxWidgets via
[MacPorts](https://www.macports.org/) as we have had some difficulty
building a wxWidgets that will work with `gmtkViz`. Note that a verfsion of GMTK is
available through MacPorts as well, so installing GMTK via MacPorts
will pull in the proper `wxWidgets`.

You will also need to install X11 on OSX 10.8 and later in order for gmtkViz
to build. We recommend installing [XQuartz](http://xquartz.macosforge.org).

`gmtkViz` is not supported on Cygwin. (If you've had any luck building
with wxWidgets on Cygwin, shoot us an email and let us in on the
secret!)


Note that the main development compiler used for GMTK is GNU
gcc/g++. Other compilers have also been tested (such as the Clang,
LLVM-GCC, and the Intel C++ compiler) but not nearly as extensively
as GNU gcc/g++. In the below, we will assume gcc/g++.

GMTK should compile with almost no warnings with the development
compilers. There are a couple of warnings in GMTK_FNGramCPT.cc that
we are looking into. 

There are a few variables that control the flags passed to the C++
compiler. You may want to specify them to the make command, though
they all have reasonable defaults. For example, you might try

```
make DEBUGFLAGS=-ggdb3 OPTFLAGS="-O3 -march=pentium4 -mfpmath=sse -ffast-math"
```

The default is just generic optimization, but you might want to use
something more aggressive (which can lead to significant speedups, but
also are such that if you run on the wrong architecture will produce
illegal instructions and/or a bus error). If you wish to turn off all
optimization, then do:

```
make OPTFLAGS= XOPTFLAGS=
```

Turning off all optimizations will lead to significantly slower
executation, but is very useful when debugging the code.


You can also change the version of gcc/g++ that you use to compile 
by specifying CC/CXX to the configure command:

```
./configure CC=gcc-9.3 CXX=g++-9.3
```

You can set the CC and CXX variables in the make command as well, but
it's safer to do so at configure time.


If you wish to create static binaries (i.e., ones that are not dependent
on shared libraries), then assuming you've got the static libraries
installed, do

```
./configure --disable-gmtkViz && make LDFLAGS=-static
```

Note, many of the above options can be combined.

Don't forget to read the [documentation](./documentation.pdf) on how
to use GMTK and on dynamic graphical models.


-- Jeff Bilmes

For general questions:  gmtk-users@u.washington.edu
Subscribe to the email list at
```
http://mailman.u.washington.edu/mailman/listinfo/gmtk-users
```
