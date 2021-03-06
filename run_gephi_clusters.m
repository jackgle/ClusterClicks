function [nodeAssign,excluded,rankAssign] = run_gephi_clusters(tempN,...
    minClust,modularity,pgRnkPrctile,javaPathVar,classPathVar,...
    toolkitPath,gv_file)

filename = 'io_test4.csv';
[status,result] = system(sprintf('"%s" -Xms512m -Xmx6144m -Dfile.encoding=Cp1252 -classpath %s;%s ClusterGephi.ClusterMain "%s" "%g" "%g"',...
    javaPathVar,classPathVar,toolkitPath,gv_file,modularity,pgRnkPrctile/100));
if status~=0
    error(result)
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