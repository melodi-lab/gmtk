function pfid = pfopen(fname, mode)
%pfid = pfopen(fname, mode) opens fname pfile for mode access.  
%returns a pfile fid, which is used with pfread/pfwrite/pfopen/pfclose
%
%This is only partially implemented and flag must be 'w' meaning write.
%mode=='w':
%If fname exists, it is overwritten. 
%The file MUST be closed with pfclose to actually write the pfile header.
%
% 2010Mar15 Arthur Kantor
if mode ~= 'w' && mode ~= 'r'
    error(['mode=' mode ' is not implemented']);
end
pfid.mode=mode;
pfid.fid=fopen(fname,mode,'b');

%I don't know the full pfile spec, so I will just copy the format
%output by feacat as closely as possible
%Also, the read mode uses the original pfinfo from Dan Ellis, and probably
%could be cleaned up 

% Basic size of pfile data items
pfid.wdsize = 4;	% bytes
pfid.version =0;

if mode == 'w'
    pfid.hdrsize = 32768;	% default, 
    pfid.num_sentences=0;
    pfid.sentStarts=zeros(64,1); %a list of 0-based row indexes for each utterance
                                 %pfid.sentStarts[pfid.num_sentences+1] refers
                                 %the next utterance to write
    pfid.num_frames=0;
    %The remaining fields get filled in on first pfwrite call

    % zero out until the start of data
    fwrite(pfid.fid, zeros(pfid.hdrsize,1), 'uchar');
    
elseif mode == 'r'
    [utrowsAll,sents,ftrs,labs,absrowAll,hdrSize] = pfinfo(fname,-1);
    pfid.hdrsize = hdrSize;
    pfid.num_sentences=sents;
    pfid.sentStarts=absrowAll;
    pfid.sentStarts(end+1)=pfid.sentStarts(end)+utrowsAll(end);
    pfid.num_frames=sum(utrowsAll);
    pfid.first_feature_column=2;
    pfid.num_features = ftrs;
    pfid.first_label_column = ftrs+2;
    pfid.num_labels = labs;
    %pfid.format = 
    
end

%save the structure
global pfidHash;
if isempty(pfidHash)
    pfidHash={};
end

pfidHash{end+1}=pfid;
pfid=numel(pfidHash);
%now we are ready to write data
