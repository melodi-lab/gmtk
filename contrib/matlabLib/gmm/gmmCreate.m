function gmm=gmmCreate(pq,q,w,m,v)
% gmm=gmmCreate(pq,q,w,m,v) create a structure representing p(q) and gaussian mixture model p(o|q)
% which is used by gmmLikelihood() to evaluate p(o,q)
% pq Qx1 is the discrete pdf p(q), where q are the states 1,2,...,Q
% q,w,m,s specify a gaussian mixture p(o|q)
% q,w,m and s all have M cols, with each ith column denoting a single component of
% a gaussian mixture, for a single discrete state q
% q is an 1xM vector over values 1,2,...Q indicating the state to
% which the ith component belongs: the component is part of p(o,q(i))
% w is an 1xM vector where w is the weight of a component for the 
% It is required that sum(w(q==j))==1 for all j in 1...Q
% m and v are DxM matrices, with ith col specifing the mean and the
% diagonal covariance of ith component.
%
%Note that you can cluster q states together by changing only pq and q.
%
%gmmCreate precomputes certain values to make calculation of p(o,q) fast

%from tuning on a intel core2 duo
blockSize=12000;


[D M]=size(m);
assert(all([D M]==size(v)));
assert(all([1 M]==size(w)));
assert(all([1 M]==size(q)));

Q = numel(unique(q));
assert(all([Q 1] == size(pq)));
assert(approxeq(sum(pq),1,1e-14));

%ginfo();

gmm.q=bsxfun(@eq,sparse(1:Q).',q ); %matrix multiplication by a sparse logical array is really fast
assert(logical(approxeq(gmm.q*w',ones(Q,1),1e-4)));

gmm.pq=pq;

gmm.maxBlock=round(blockSize/D);
gmm.maxBlock=min(gmm.maxBlock,M);

gmm.m=[m m(:,1:(gmm.maxBlock-1))];
gmm.invVar=-.5./[v v(:,1:(gmm.maxBlock-1))];

gmm.addCoef= log(pq(q,:))+log(w')-log(2*pi)*D/2-sum(log(v))'/2;
