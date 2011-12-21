function add_diagonal_covariance_component_in_full_cov_form(GMTK_means_filename, GMTK_diag_covs_filename, GMTK_dlinks_filename, ngc, ndim, laa, diag_cov_weight)

if (isa(ngc, 'char'))
        ngc = str2num(ngc);
end

if (isa(ndim, 'char'))
        ndim = str2num(ndim);
end

if (isa(laa, 'char'))
        laa = str2num(laa);
end

if (isa(diag_cov_weight, 'char'))
        diag_cov_weight = str2num(diag_cov_weight);
end

format long;

fid_means = fopen(GMTK_means_filename, 'r');

if (fid_means == -1)

	'Error: Could not open the GMTKs learned mean vectors file'
	return
end

a_line = fgets(fid_means);
a_line = fgets(fid_means);
number_of_gc=fscanf(fid_means, '%d');

GMTK_means = zeros(ndim, ngc);

for j = 1:ngc

	index_value=fscanf(fid_means, '%d');
	a_line = fgets(fid_means);
	[token, a_line] = strtok(a_line);
	[token, a_line] = strtok(a_line);
	
	for k = 1:ndim

		[token, a_line] = strtok(a_line);
		GMTK_means(k, j) = str2num(token); 
	end	
end

fclose(fid_means);

fid_covs = fopen(GMTK_diag_covs_filename, 'r');

if (fid_covs == -1)

        'Error: Could not open the GMTKs learned diagonal covariance vectors file'
        return
end

a_line = fgets(fid_means);
a_line = fgets(fid_means);
number_of_gc=fscanf(fid_covs, '%d');

GMTK_diag_covs = zeros(ndim, ngc);

for j = 1:ngc

        index_value=fscanf(fid_covs, '%d');
        a_line = fgets(fid_covs);
        [token, a_line] = strtok(a_line);
        [token, a_line] = strtok(a_line);

        for k = 1:ndim

                [token, a_line] = strtok(a_line);
                GMTK_diag_covs(k, j) = str2num(token);
        end
end

fclose(fid_covs);

fid_dlinks = fopen(GMTK_dlinks_filename, 'r');
 
if (fid_dlinks == -1) 
 
        'Error: Could not open the GMTKs learned dlinks file'
        return 
end 
 
a_line = fgets(fid_dlinks);
a_line = fgets(fid_dlinks);

number_of_gc=fscanf(fid_dlinks, '%d');

GMTK_dlinks = zeros(ndim, (laa+1)*ndim, ngc); 

for j = 1:ngc

    index_value=fscanf(fid_dlinks, '%d');
    dlink_mat_str = fgets(fid_dlinks);
    dlink_structure_str = fgets(fid_dlinks);
    ndim_str = fgets(fid_dlinks);

    for k = 1:ndim

        a_line = fgets(fid_dlinks); 
	    [token, a_line] = strtok(a_line); 
     	ndlinks = str2num(token); 
		
        for m = (k+1):(ndlinks+k)

            [token, a_line] = strtok(a_line);
      	    GMTK_dlinks(k, m, j) = str2num(token);
        end
    end
end

fclose(fid_dlinks);

%GMTK_means
%GMTK_diag_covs
%GMTK_dlinks(3,18,:)

Bs = GMTK_dlinks; %not a square matrix in general
B1s = zeros(ndim, ndim, ngc);
B2s = zeros(ndim, laa*ndim, ngc);

U1s = zeros(ndim, ndim, ngc);
U2s = zeros(ndim, ndim*laa, ngc);
Ds = zeros(ndim, ndim, ngc);

for j = 1:ngc
	
	B1s(:,:,j) = Bs(:,1:ndim,j);
	B2s(:,:,j) = Bs(:,ndim+1:end, j);
	Ds(:,:,j) = inv(diag(GMTK_diag_covs(:,j)));
	U1s(:,:,j) = eye(ndim) - B1s(:,:,j);
	U2s(:,:,j) = inv(U1s(:,:,j))*B2s(:,:,j);
end

means = zeros(ndim, ngc);
covs = zeros(ndim, ndim, ngc);
Ks_new = zeros(ndim, ndim, ngc);
Rs_new = zeros(ndim, ndim, ngc);
D12s_new = zeros(ndim, ndim, ngc);
U1s_new = zeros(ndim, ndim, ngc);
B1s_new = zeros(ndim, ndim, ngc);
B2s_new = zeros(ndim, laa*ndim, ngc);
Ds_new = zeros(ndim, ndim, ngc);
GMTK_means_new = zeros(ndim, ngc);
GMTK_diag_covs_new = zeros(ndim, ngc);
GMTK_dlinks_new = zeros(ndim, (laa+1)*ndim, ngc);

for j = 1:ngc

	means(:,j) = inv(U1s(:,:,j))*GMTK_means(:,j);
	covs(:,:,j) = inv(U1s(:,:,j)'*Ds(:,:,j)*U1s(:,:,j));
	covs(:,:,j) = (1-diag_cov_weight)*covs(:,:,j) + diag_cov_weight*eye(ndim); %add a diagonal covariance component here
	Ks_new(:,:,j) = inv(covs(:,:,j));
	Rs_new(:,:,j) = chol(Ks_new(:,:,j));
	D12s_new(:,:,j) = diag(diag(Rs_new(:,:,j)));
	U1s_new(:,:,j) = D12s_new(:,:,j)\Rs_new(:,:,j);
	Ds_new(:,:,j) = D12s_new(:,:,j).^2;	
	GMTK_means_new(:,j) = U1s_new(:,:,j)*means(:,j);
	GMTK_diag_covs_new(:,j) = diag(inv(Ds_new(:,:,j)));
	B1s_new(:,:,j) = eye(ndim) - U1s_new(:,:,j);
	B2s_new(:,:,j) = U1s_new(:,:,j)*U2s(:,:,j);
	GMTK_dlinks_new(:,:,j) = [B1s_new(:,:,j),B2s_new(:,:,j)] ;	
end 

fid = fopen(GMTK_means_filename, 'w');

fprintf(fid, '\n');
fprintf(fid, '%d\n', ngc);

for j = 1:ngc

	fprintf(fid, '%d\n', (j-1));
	fprintf(fid, '%s%s ', 'meanEC', num2str(j-1));
	fprintf(fid, '%d ', ndim);

	for k=1:ndim
        	fprintf(fid, '%.11f ', GMTK_means_new(k,j));
	end

	fprintf(fid, '\n');
end

fclose(fid);

fid = fopen(GMTK_diag_covs_filename, 'w');

fprintf(fid, '\n');
fprintf(fid, '%d\n', ngc);

for j = 1:ngc

	fprintf(fid, '%d\n', (j-1));
	fprintf(fid, '%s%s ', 'covarEC', num2str(j-1));
	fprintf(fid, '%d ', ndim);

	for k=1:ndim

        	fprintf(fid, '%.15f ', GMTK_diag_covs_new(k,j));
	end

	fprintf(fid, '\n');
end

fclose(fid);

fid = fopen(GMTK_dlinks_filename, 'w');

fprintf(fid, '\n');
fprintf(fid, '%d\n', ngc);

for j = 1:ngc

	fprintf(fid, '%d\n', (j-1));
	fprintf(fid, '%s%s\n', 'dlink_mat', num2str(j-1));
	fprintf(fid, '%s%s\n', 'dlink_str', num2str(j-1));
	fprintf(fid, '%d\n', ndim);

	for k=1:ndim

        	fprintf(fid, '%d ', ((laa+1)*ndim-k));

        	for m = (k+1):(laa+1)*ndim

                	fprintf(fid, '%.15f ', GMTK_dlinks_new(k,m,j));
        	end

        	fprintf(fid, '\n');
	end
end

