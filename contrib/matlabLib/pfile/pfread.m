function [D,L,I] = pfread(N,U,S,F)
% [D,L,I] = pfread(N,U,S,R)  Read one utterance from a pfile
%    Read the data and labels for utterance U from pfile named N.
%    S is starting frame (optional, default 0) and F is the 
%    maximum number of frames to read (optional, default to EOS). 
%    Data is returned in D, labels in L.
%    I is a nx2 matrix containing the [uttId frameId] for each returned frame
%
%    N can also be a pfid opened with pfopen()
%    this is a lot faster for repeated reads from the same file
%
%    U can also be a 1x2 matrix specifying a range of utterances which are
%    read and returned in one call
%    
% 1998aug14 dpwe@icsi.berkeley.edu   M-script version
% 2007sept16 akantor@uiuc.edu can also use pfid from pfopen() instead of
% filename

if nargin < 3
  S = 0;
end
if nargin < 4
  F = -1;
end

if numel(U)==1
    U=[U U];
else
    assert(nargin < 3,'cannot read multiple partial utterances')
end
U=U+1;

global pfidHash;
if ischar(N)
    % Use our sister function to get the raw index into the pfile, from 
    % the header and index
    
    [utrowsAll,dummy,ftrs,labs,absrowAll] = pfinfo(N,-1);
    
    fid = fopen(N,'r','b');
    utrows=sum(utrowsAll(U(1):U(2)));
    absrow=absrowAll(U(1));
else
    %persistent pfid utrowsAll sents ftrs labs absrowAll ;
    assert(pfidHash{N}.mode =='r', 'pfid is not opened in read mode')
    fid=pfidHash{N}.fid;
    utrows=pfidHash{N}.sentStarts(U(2)+1)-pfidHash{N}.sentStarts(U(1));
    absrow=pfidHash{N}.sentStarts(U(1));
    ftrs=pfidHash{N}.num_features;
    labs=pfidHash{N}.num_labels;
end


%disp('opened big endian');
% Seek to the row
% Assume std header size
hdrsize=32768;
wdsize=4;
prfxcols = 2;	% 2 columns of prefix (utt, frm) on each pfile row
rowlen = prfxcols+ftrs+labs;

D = [];
L = [];

% Which rows are we actually going to read?
startrow = min(S, utrows);
nrows = utrows - startrow;
if F >= 0
  nrows = min(F, nrows);
end

% First read feature data, if any
if nargout>2
    fseek(fid, hdrsize+wdsize*rowlen*(absrow+startrow), -1);
    I = fread(fid, [prfxcols,nrows], [num2str(prfxcols) '*int32'], (ftrs+labs)*wdsize)';
end

if ftrs > 0
	fseek(fid, hdrsize+wdsize*rowlen*(absrow+startrow), -1);
	D = fread(fid, [prfxcols+ftrs+labs,nrows], 'float32');
	% collapse away the non-feature cols (& transpose)
	D = D(prfxcols+(1:ftrs),:)';
end

% Read label data, if any
if labs > 0
	fseek(fid, hdrsize+wdsize*rowlen*(absrow+startrow), -1);
	L = fread(fid, [prfxcols+ftrs+labs,nrows], 'int32');
	% collapse away the non-label cols (& transpose)
	L = L(prfxcols+ftrs+(1:labs),:)';
end

% done
if ischar(N)
    fclose(fid);
end
