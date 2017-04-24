[specNorm,diffNormSpec] = spec_norm_diff(sumSpec,p.stIdx,p.edIdx);

excludedIn = 1:length(dTTmatNorm);

[specDist,~,~] = spectra_dist(diffNormSpec(excludedIn,p.stIdx:p.edIdx));


[iciDist,~,~,~] = ici_dist_mode(dTTmatNorm(excludedIn,minICI:maxICI),...
    p.barInt(minICI:maxICI));
compDist = squareform(specDist.*iciDist,'tomatrix');
wScore = [];
bScore = [];
for iN = 1:length(naiItr)
    wInDist = [];
    btwnDist = [];
    for iK = 1:length(naiItr{iN})
        notThisSet = setdiff(excludedIn,naiItr{iN}{iK});
        wInDist(iK)= sum(sum(compDist(naiItr{iN}{iK},naiItr{iN}{iK}).^2))/2;
        btwnDist(iK) =  sum(sum(compDist(naiItr{iN}{iK},notThisSet).^2))/2;
    end
    wScore(iN,1) = sum(wInDist);
    bScore(iN,1) = sum(btwnDist);
end
dispMat = [ka,wScore./10000,bScore./10000,(wScore-bScore)./10000]