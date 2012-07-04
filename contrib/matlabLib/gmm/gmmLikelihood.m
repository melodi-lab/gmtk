function [like,maxmah]=gmmLikelihood(o,gmm)
% [like,maxmah]=gmmLikelihood(o,gmm) calculate the likelihood p(o,q) using the mixture of gaussians model
% o is an DxN matrix of N D-dimentional observations
% gmm is the mixture of gaussians model created with gmmCreate(). gmm
% components should work with D-dimentional observations
% output is an QXN matrix where Q is cardinality of q
% 
% see gmmCreate()

% ### uncomment this to run things in parallel
%o=distributed(o);
%spmd
%o=getLocalPart(o);

N=size(o,2);
[Q,M]=size(gmm.q);
like=zeros(Q,N);
maxmah=zeros(1,N);


%Straightforward implementation:
%process each o seperately against all components at once
% for k =1:N
%     mah=sum(bsxfun(@minus, o(:,k), gmm.m).^2.*gmm.invVar); %mahalanobis
%     l(:,k) = gmm.q*(exp(mah+gmm.addCoef).');
%     [~,idx]=max(mah);
%     maxmah(k)=idx;
% end

%an implementation to improve cpu caching: about 2.5 times faster

blockRanges=1:gmm.maxBlock:N+1;
if (blockRanges(end) ~= N+1)
    blockRanges(end+1)=N+1;
end


%process in blocks that fit in CPU cache
for b=1:numel(blockRanges)-1
    cbr=[blockRanges(b),blockRanges(b+1)-1];
    block=o(:,cbr(1):cbr(2));
    blockSz=diff(cbr)+1;
    mah=zeros(M,blockSz);
    for k =1:M
        mah(k,:)=sum((block-gmm.m(:,k:k+blockSz-1)).^2.*gmm.invVar(:,k:k+blockSz-1)); %mahalanobis
    end
    %now circular shift down the mah for each o into the right position
    for l=2:blockSz
        mah(:,l)=[mah(end-l+2:end,l);mah(1:end-l+1,l)];
    end
    [dummy,idx]=max(mah(:,1:blockSz));
    maxmah(cbr(1):cbr(2))=idx;
    %FIXME store liklihoods as log-liklihoods to preserve precision in
    %32-bit floats
    like(:,cbr(1):cbr(2)) = gmm.q*(exp(bsxfun(@plus,mah(:,1:blockSz),gmm.addCoef)));
end
% ### uncomment this to run things in parallel
%end
%like=[like{:}];
%maxmah=[maxmah{:}];


