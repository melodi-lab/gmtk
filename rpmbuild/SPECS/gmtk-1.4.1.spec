%define debug_package %{nil}

Name:		gmtk
Version:	1.4.1
Release:	1%{?dist}
Summary:	The Graphical Models Toolkit
BuildRoot:	/var/tmp/rpmbuild/gmtk-1.4.1-root
Group:		Applications/Engineering
License:	OSL-3.0
URL:		https://j.ee.washington.edu/trac/gmtk
Source:		http://melodi.ee.washington.edu/downloads/gmtk/gmtk-1.4.1.tar.gz

AutoReq: no
BuildRequires:	hdf5-devel wxGTK3-devel libstdc++-devel
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
%configure --with-wx-config=wx-config-3.0
make %{?_smp_mflags}

%install
%makeinstall

%files
%{_bindir}/discrete-mi
%{_bindir}/fixTri.sh
%{_bindir}/generate_random_graph.pl
%{_bindir}/gmtkDMLPtrain
%{_bindir}/gmtkDTindex
%{_bindir}/gmtkEMtrain
%{_bindir}/gmtkJT
%{_bindir}/gmtkKernel
%{_bindir}/gmtkModelInfo
%{_bindir}/gmtkMMItrain
%{_bindir}/gmtkNGramIndex
%{_bindir}/gmtkOnline
%{_bindir}/gmtkParmConvert
%{_bindir}/gmtkPrint
%{_bindir}/gmtkTFmerge
%{_bindir}/gmtkTie
%{_bindir}/gmtkTime
%{_bindir}/gmtkTriangulate
%{_bindir}/gmtkViterbi
%{_bindir}/gmtkViz
%{_bindir}/obs-cat
%{_bindir}/obs-concat
%{_bindir}/obs-diff
%{_bindir}/obs-info
%{_bindir}/obs-print
%{_bindir}/obs-skmeans
%{_bindir}/obs-stats
%{_bindir}/obs-window
%{_bindir}/triangulateGA
%{_bindir}/triangulateParallel
%{_bindir}/triangulateTimings
%doc COPYING
%doc NEWS
%{_mandir}/man1/*

%changelog
* Fri Oct 02 2015 Richard Rogers <rprogers@uw.edu> 1.4.1-1
- Simplified specfile
- Added man pages
