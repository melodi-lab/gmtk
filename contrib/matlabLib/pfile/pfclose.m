function pfclose(pfid)
% pfclose(pfid) closes an open pfile.
%This is especially necessary when writing to pfile because the header is
%written by this function

global pfidHash;

fid=pfidHash{pfid}.fid;
if pfidHash{pfid}.mode=='w'
    
    assert(pfidHash{pfid}.num_sentences>0,'No utterances written: cannot save an empty pfile.');
    % write the row index for each utterance, plus the start of the next
    % yet-unwritten utterance
    fwrite(fid,pfidHash{pfid}.sentStarts(1:pfidHash{pfid}.num_sentences+1), 'int32',0);
    
    %write the header
    num_frames=pfidHash{pfid}.sentStarts(pfidHash{pfid}.num_sentences+1);
    fseek(fid, 0, -1);
    ncol=2+pfidHash{pfid}.num_features+pfidHash{pfid}.num_labels;
    fprintf(fid,'-pfile_header version %d size %d\n',pfidHash{pfid}.version, pfidHash{pfid}.hdrsize);
    fprintf(fid,'-num_sentences %d\n',pfidHash{pfid}.num_sentences);
    fprintf(fid,'-num_frames %d\n',num_frames);
    fprintf(fid,'-first_feature_column %d\n',pfidHash{pfid}.first_feature_column);
    fprintf(fid,'-num_features %d\n',pfidHash{pfid}.num_features);
    fprintf(fid,'-first_label_column %d\n',pfidHash{pfid}.first_label_column);
    fprintf(fid,'-num_labels %d\n',pfidHash{pfid}.num_labels);
    fprintf(fid,'-format %s\n',pfidHash{pfid}.format);
    fprintf(fid,'-data size %ld offset 0 ndim 2 nrow %d ncol %d\n', num_frames*ncol, num_frames, ncol);
    fprintf(fid,'-sent_table_data size %d offset %ld ndim 1\n',pfidHash{pfid}.num_sentences+1, num_frames*ncol);
    fprintf(fid,'-end\n');
end

fclose(fid);

pfidHash{pfid}=[];