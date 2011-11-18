function pfwrite(pfid, o, id)
%pfwrite(pfid, o) writes o to pfile handle pfid
%o is a DxM matrix, with M frames and each frame with D float features.
%All pfwrites to pfid except the first must have o with the same number of
%features D, otherwise an error is raised.
%
% If id is specified, it is a 2xM matrix, listing the sentenceId and frame
% Id in the first and second rows.  It can be used to write multiple
% sentences in one pfwrite call. The utterance ids must be sequential.
% They are adjusted to be consecutive with the previously written
% utterances.
% The frame ids must be 0-based and sequential within each utterance.
% The validity of the frame ids is not checked.
%
% If id is not give, o is assumed to be a single utterance.
%
% writing labels is not yet implemented, mainly because other tools
% (feacat, obs-print) have problems with pfiles containing features and
% labels
%
% 2010Mar15 Arthur Kantor

global pfidHash;

[D,N]=size(o);
if pfidHash{pfid}.num_sentences>0
    assert(pfidHash{pfid}.num_features == D, ...
        sprintf('number of features is inconsistent. Should be %d',pfidHash{pfid}.num_features));
else
    pfidHash{pfid}.first_feature_column=2;
    pfidHash{pfid}.num_features = D;
    pfidHash{pfid}.first_label_column = D+2;
    pfidHash{pfid}.num_labels = 0;
    pfidHash{pfid}.format = ['dd' repmat('f',1,D)];
end

%prefix the utt sent to each frame

if nargin<3
    id=[repmat(pfidHash{pfid}.num_sentences,1,N); 0:(N-1)];
    newStarts=N;
else
    d=diff(id(1,:));
    d(end+1)=1;
    assert(all(d==0 | d==1),'sentence Ids must be consecutive');
    newStarts=find(d);
    id(1,:)=id(1,:)-id(1,1)+pfidHash{pfid}.num_sentences;
end
addedSents=numel(newStarts);
if pfidHash{pfid}.num_sentences+addedSents > numel(pfidHash{pfid}.sentStarts)
    pfidHash{pfid}.sentStarts=[pfidHash{pfid}.sentStarts;zeros(pfidHash{pfid}.num_sentences+addedSents,1)];
end
curEnd=pfidHash{pfid}.sentStarts(pfidHash{pfid}.num_sentences+1);
pfidHash{pfid}.sentStarts(pfidHash{pfid}.num_sentences+1+(1:addedSents))=curEnd+newStarts;
pfidHash{pfid}.num_sentences=pfidHash{pfid}.num_sentences+addedSents;


%now write the actual utterances 
% for i=1:N
%     fwrite(pfidHash{pfid}.fid, id(:,i), 'int32',0);
%     fwrite(pfidHash{pfid}.fid, o(:,i), 'float32',0);
% end
% The following is a faster implementation
fwrite(pfidHash{pfid}.fid, [reshape(typecast(int32(id(:)),'single'),2,[]);cast(o,'single')], 'single',0);
