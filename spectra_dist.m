function [distClickE,rows1,cols1,specClickTf_norm,specClickTf_norm_short] = ...
    spectra_dist(specClickTf,stIdx,edIdx)

% TODO: make weighted correlation an option to be fed in
wcorTF = 0;

minSSsection = min(specClickTf,[],2);
specClickTf_minNorm = (specClickTf - ...
  minSSsection(:,ones(1,size(specClickTf,2))));
maxSSsection = max(specClickTf_minNorm,[],2);
specClickTf_norm = specClickTf_minNorm./maxSSsection(:,ones(1,size(specClickTf_minNorm,2)));
specClickTf_norm_short = specClickTf_norm(:,stIdx:edIdx);
specClickTf_norm_short = diff(specClickTf_norm_short,1,2);

tempN = size(specClickTf_norm_short,1);
offaxN = ((tempN.^2)-tempN)./2;

rows1 = zeros(offaxN, 1);
cols1 = zeros(offaxN, 1);
n = 1;

for itrA = 1:size(specClickTf_norm_short,1)-1  
    for itrB = itrA+1:size(specClickTf_norm_short,1)
        rows1(n) = itrA;
        cols1(n) = itrB;
        n = n+1;
    end
end
if ~wcorTF
    % distClick = pdist(specClickTf_norm,'seuclidean');
    distClick = pdist(specClickTf_norm_short,'correlation');
    distClickE = (exp(-distClick));
else
    % weighted correlation
    % coordinate weights    
    %Wgts = 1:-1/(size(specClickTf,2)-1):0; 
    Wgts = ones(1,size(specClickTf,2));   
    % coordinate weights
    X = specClickTf_norm;
    X = bsxfun(@minus,X,mean(X,2));
    Xmax = max(abs(X),[],2);
    X2 = bsxfun(@rdivide,X,Xmax);
    Xnorm = sqrt(sum(X2.^2, 2));
    Xnorm(Xmax==0) = 0;
    Xnorm = Xnorm .* Xmax;
    X = bsxfun(@rdivide,X,Xnorm);
    
    % weuc = @(XI,XJ,W)(sqrt(bsxfun(@minus,XI,XJ).^2 * W'));
    wcorr= @(XI,XJ,W)(1-sum(bsxfun(@times,XI,XJ) * W',2));
    Dwgt = pdist(X, @(Xi,Xj,Wgts) wcorr(Xi,Xj,Wgts),Wgts);
    
    distClickE = (exp(-real(Dwgt)));
end
