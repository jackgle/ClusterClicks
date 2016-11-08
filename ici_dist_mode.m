function [distClickE,rows1,cols1,iciMode] = ici_dist_mode(iciMat,barInt,minBinIdx)
% iciNorm = iciMat./repmat(sum(iciMat,2),1,size(iciMat,2));

[~,iciModeIdx] = max(iciMat(:,minBinIdx:end),[],2);
iciMode = barInt(iciModeIdx+minBinIdx-1) + barInt(2)./2;

tempN = size(iciMode,1);
offaxN = ((tempN.^2)-tempN)./2;

rows1 = zeros(offaxN, 1);
cols1 = zeros(offaxN, 1);
n = 1;

for itrA = 1:size(iciMode,1)-1  
    for itrB = itrA+1:size(iciMode,1)
        rows1(n) = itrA;
        cols1(n) = itrB;
        n = n+1;
    end
end
distClick = pdist(iciMode','euclidean');

distClickE = exp(-distClick);

