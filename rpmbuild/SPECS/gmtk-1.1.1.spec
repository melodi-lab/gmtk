
Name:		gmtk
Version:	1.1.1
Release:	1%{?dist}
Summary:	The Graphical Models Toolkit
BuildRoot:	/var/tmp/rpmbuild/gmtk-1.1.1-root
Group:		Applications/Engineering
License:	OSL-3.0
URL:		https://j.ee.washington.edu/trac/gmtk
Source:		http://melodi.ee.washington.edu/downloads/gmtk/gmtk-1.1.1.tar.gz

BuildRequires:	hdf5-devel wxGTK3-devel
Requires:	hdf5 wxGTK3

%description
The Graphical Models Toolkit (GMTK) is an open source, publicly available 
toolkit for rapidly prototyping statistical models using dynamic graphical 
models (DGMs) and dynamic Bayesian networks (DBNs). GMTK can be used for 
applications and research in speech and language processing, bioinformatics, 
activity recognition, and any time series application. GMTK has many features, 
including exact and approximate inference; a large variety of built-in factors 
including dense, sparse, and deterministic conditional probability tables, 
native support for ARPA backoff-based factors and factored language models, 
parameter sharing, gamma and beta distributions, dense and sparse Gaussian 
factors, heterogeneous mixtures, deep neural network factors, and 
time-inhomogeneous trellis factors; arbitrary order embedded Markov chains; 
a GUI-based graph viewer; flexible feature-file support and processing tools 
(supporting pfiles, HTK files, ASCII/binary, and HDF5 files); and both offline 
and streaming online inference methods that can be used for both parameter 
learning and prediction. More information is available in the documentation. 
All in all, GMTK offers a flexible, concise, and expressive probabilistic 
modeling framework with which one may rapidly specify a vast collection of 
temporal statistical models.

%prep
%setup -q


%build
%configure
make %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT


%clean
rm -rf $RPM_BUILD_ROOT


%files
%defattr(-,root,root,-)
/usr/bin/discrete-mi
/usr/bin/fixTri.sh
/usr/bin/generate_random_graph.pl
/usr/bin/gmtkDMLPtrain
/usr/bin/gmtkDTindex
/usr/bin/gmtkEMtrain
/usr/bin/gmtkJT
/usr/bin/gmtkKernel
/usr/bin/gmtkModelInfo
/usr/bin/gmtkNGramIndex
/usr/bin/gmtkOnline
/usr/bin/gmtkParmConvert
/usr/bin/gmtkPrint
/usr/bin/gmtkTFmerge
/usr/bin/gmtkTie
/usr/bin/gmtkTime
/usr/bin/gmtkTriangulate
/usr/bin/gmtkViterbi
/usr/bin/gmtkViz
/usr/bin/obs-cat
/usr/bin/obs-concat
/usr/bin/obs-diff
/usr/bin/obs-info
/usr/bin/obs-print
/usr/bin/obs-skmeans
/usr/bin/obs-stats
/usr/bin/triangulateGA
/usr/bin/triangulateParallel
/usr/bin/triangulateTimings
%doc COPYING
%doc NEWS

%changelog
