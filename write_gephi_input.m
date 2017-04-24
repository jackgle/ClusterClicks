function [file2Write,isolated] = write_gephi_input(distMat, minClust,thr)
excluded = [];

tempN = size(distMat,1);
offaxN = ((tempN.^2)-tempN)./2;
rows1 = zeros(offaxN, 1);
cols1 = zeros(offaxN, 1);
distMatL = zeros(offaxN, 1);
n = 1;

for itrA = 1:tempN-1
    for itrB = itrA+1:tempN
        rows1(n) = itrA;
        cols1(n) = itrB;
        distMatL(n) = distMat(itrA,itrB);
        n = n+1;
    end
end

if isempty(thr)
    thrP = prctile(distMatL,90);
else
    thrP = prctile(distMatL,thr);
end
tempN = ceil(sqrt(length(distMatL)*2));
% [values, rows, cols, isolated] = prune_by_node1(tempN, distMatL, rows1, cols1, thrP);
[values, rows, cols, isolated] = prune(tempN, distMatL, rows1, cols1, thrP);
nodesRemaining = length(values);
if nodesRemaining < minClust
    fprintf('Too few nodes remaining after pruning, skipping iteration \n')
    excluded = isolated;
    return % if there aren't nodes/edges left, skip iteration
end

file2Write = sprintf('click_output.gv');
fid = fopen(file2Write, 'W');
fprintf(fid, 'strict graph {\n');
for vidx=1:length(values)
    fprintf(fid, 'n%d -- n%d [weight=%f, style=invis];\n', ...
        rows(vidx),cols(vidx), values(vidx));
end
fprintf(fid, '}\n');
fclose(fid);
