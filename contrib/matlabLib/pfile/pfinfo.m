function [R,U,F,L,T,hdrsize] = pfinfo(N,O)
% [R,U,F,L,T] = pfinfo(name[,utt])   Return info for a pfile
%    Open an ICSI feature file ("pfile") named <name>.  Return the 
%    number of utterances in U, the features per frame in F, the labels 
%    per frame in L and the total number of rows in R.  If <utt> is 
%    specified, R is the number of rows in that utterance only, and T
%    is the row number where it starts (zero otherwise).
% if O == -1, the R and T are vectors for all U utterances
% 1998aug14 dpwe@icsi.berkeley.edu  M-script version
% 2007sept16 akantor@uiuc.edu vectorized R and T and returns the header size

if nargin < 2
	O = -2;
end

fid = fopen(N,'r','b');
if fid < 0
	error(['couldnt open pfile ', N]);
end

% Basic size of pfile data items
wdsize = 4;	% bytes

% Read through the pfile header
done = 0;			% flag to end loop when -end read
first = 1;			% semaphore to check header starts as expected
ixoffs = -1;		% flag that no sent_table_data found
hdrsize = 32768;	% default, but will read val from header
while done == 0 
	s = fscanf(fid, '%s', 1);
	if first == 1
		if strcmp(s, '-pfile_header') == 0
			fclose(fid)
			error([N,' is not a pfile (', s,')'])
		end
		block = '-pfile_header';
		first = 0;
	elseif strcmp(s, '-end') == 1
		done = 1;
	elseif strcmp(s, '-data')
		block = '-data';
	elseif strcmp(s, '-sent_table_data')
		block = '-sent_table_data';
        % put other argless tags in here
	else
		% Assume it is followed by a val
		val = fscanf(fid, '%s', 1);
		% Maybe pick it out now
		if strcmp(s, '-num_sentences') == 1
			U = str2num(val);
		elseif strcmp(s, '-num_frames') == 1
			R = str2num(val);
		elseif strcmp(s, '-num_features') == 1
			F = str2num(val);
		elseif strcmp(s, '-num_labels') == 1
			L = str2num(val);
		elseif strcmp(block, '-sent_table_data') == 1 & strcmp(s, 'offset') == 1
			ixoffs = str2num(val);
		elseif strcmp(block, '-pfile_header') == 1 & strcmp(s, 'size') == 1
			hdrsize = str2num(val);
		end
	end
end

T = 0;

if O >= -1
	% Read the sent_table_data to get data on one seg
	if ixoffs < 0
		fclose(fid);
		error('specific segment info requested, but pfile has no index');
	end
	% seek to the sentence table data
	fseek(fid, hdrsize + wdsize*ixoffs, -1);
	% read the row index for each utterance
	tab = fread(fid, inf, 'int32');
	if (size(tab,1) ~= (U+1))
		fclose(fid);
		error([num2str(U), ' sentences but got ', num2str(size(tab,1)), ' index entries']);
    end
    if O >= 0
    	% return data for requested utterance
    	T = tab(1+O);
    	R = tab(1+O+1) - T;
    else
        %return data for all utterances
        T=tab(1:end-1);
        R = diff(tab);
    end
end

fclose(fid);

		
