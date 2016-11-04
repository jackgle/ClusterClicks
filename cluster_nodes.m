function [nodeAssign,excluded,rankAssign] = cluster_nodes(distMat,...
    minClust,thr,modularity,pgRnkPrctile)
excluded = [];
nodeAssign = {};

filename = 'io_test4.csv';

tempN = size(distMat,1);
offaxN = ((tempN.^2)-tempN)./2;

rows1 = zeros(offaxN, 1);
cols1 = zeros(offaxN, 1);
distMatL = zeros(offaxN, 1);
n = 1;

for itrA = 1:size(distMat,1)-1
    for itrB = itrA+1:size(distMat,1)
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

javaPathVar = 'C:\Program Files\Java\jre6\bin\java.exe';
classPathVar = ' E:\workspace\ClusterGephi_sio\bin;';
toolkitPath = 'E:\workspace\ClusterGephi_sio\gephi-toolkit-0.8.7-all\gephi-toolkit.jar';
[status,result] = system(sprintf('"%s" -Xms512m -Xmx6144m -Dfile.encoding=Cp1252 -classpath %s;%s ClusterGephi.ClusterMain "%s" "%g" "%g"',...
    javaPathVar,classPathVar,toolkitPath,file2Write,modularity,pgRnkPrctile));
if status~=0
    error(result)
% else
%     disp(result) 
end
% Read in output table from csv created by gephi
[dataTable,labels] = xlsread(filename);
% pgRank = dataTable(3,:)>=prctile(dataTable(3,:),50);
clusters = dataTable(3,:)';
nodeNums = char(labels(1,:)');
rankVals = dataTable(2,:);
nodeNums = str2double(cellstr(nodeNums(:,2:end)));
excluded = setdiff(1:tempN,nodeNums);
clustBins = 0:max(clusters);
counts = histc(clusters,clustBins);
keepClust = find(counts >= minClust);
clustNums = clustBins(keepClust);

if ~isempty(keepClust)
    for i4 = 1:length(keepClust)        
        nodeAssign{i4} = nodeNums(clusters==clustNums(i4));
        rankAssign{i4} = rankVals(clusters==clustNums(i4));
    end
end
%excluded = setdiff(1:size(distMat,1),cell2mat(nodeAssign'));
1;