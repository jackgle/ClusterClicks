function [spectraMean,clickAssign,clustSizes,spectraHolder,...
    clickAssignAll] = cluster_clicks(specClickTf,p,javaPathVar,...
    classPathVar,toolkitPath)

clickAssign = [];
spectraHolder = [];
clustSizes = 0;
filename = 'io_test4.csv';
wcorTF = 0; % flag for weighted corrlation. 
% Weighted option lets you give more weight to similarities at some
% frequencies than others.  

% Cepstral option -
% specCep = fft(specClickTrunc,[],2);
% specClickTf_nDiff = specCep(:,4:50);
% Slope normalization -
% alphaVec = specClickTrunc(:,end)-specClickTrunc(:,1);
% slopeMat = zeros(size(specClickTrunc));
% interpVal = size(specClickTrunc,2)-1;
% for i2 = 1:length(alphaVec)
%     slopeMat(i2,:) = 0:alphaVec(i2)/interpVal:alphaVec(i2);
% end
% specClickTf_SlopeNorm = (specClickTrunc - slopeMat);

% Normalize clicks on 0-1 scale
specClickTrunc = specClickTf;%(:,p.stIdx:p.edIdx);
minSSsection = min(specClickTrunc,[],2);
specClickTf_minNorm = (specClickTrunc - ...
     minSSsection(:,ones(1,size(specClickTrunc,2))));
maxSSsection = max(specClickTf_minNorm,[],2);
specClickTf_norm = specClickTf_minNorm./maxSSsection(:,ones(1,size(specClickTf_minNorm,2)));
specClickTf_norm_short = specClickTf_norm(:,p.stIdx:p.edIdx);
if p.diff
    specClickTf_norm_short = diff(specClickTf_norm_short,1,2);
end

% specDiff = diff(specClickTrunc,1,2);
% specClickTf_nDiff = specDiff./repmat(std(specDiff,[],2),1,size(specDiff,2));

% clickCepTrunc = clickCepstra(:,1:35);
% find distance between them - no warping!
tempN = size(specClickTf_norm,1);
offaxN = ((tempN.^2)-tempN)./2;

rows1 = zeros(offaxN, 1);
cols1 = zeros(offaxN, 1);
n = 1;

for itrA = 1:size(specClickTf_norm,1)-1
    for itrB = itrA+1:size(specClickTf_norm,1)
        rows1(n) = itrA;
        cols1(n) = itrB;
        n = n+1;
    end
end

if wcorTF ~= 1
    % Can use euclidean distance, but it doesn't capture shape very well
    %distClick = pdist(specClickTf_norm,'seuclidean');

    distClick = pdist(specClickTf_norm_short,'correlation');
    distClickE = (exp(-distClick));
else
    % weighted correlation
    % Wgts = 1:-1/(size(specClickTf,2)-1):0; % Vector from 1 to 0 ->
    % decreasing weight as frequency increases
    Wgts = ones(1,size(specClickTf_norm_short,2));  %all ones -> no weighting 
    
    X = specClickTf_norm_short;
    X = bsxfun(@minus,X,mean(X,2));
    Xmax = max(abs(X),[],2);
    X2 = bsxfun(@rdivide,X,Xmax);
    Xnorm = sqrt(sum(X2.^2, 2));
    Xnorm(Xmax==0) = 0;
    Xnorm = Xnorm .* Xmax;
    X = bsxfun(@rdivide,X,Xnorm);
    
    % euclidean option:
    % weuc = @(XI,XJ,W)(sqrt(bsxfun(@minus,XI,XJ).^2 * W'));
    
    % correlation option:
    wcorr= @(XI,XJ,W)(1-sum(bsxfun(@times,XI,XJ) * W',2));
    Dwgt = pdist(X, @(Xi,Xj,Wgts) wcorr(Xi,Xj,Wgts),Wgts);
    
    distClickE = (exp(-real(Dwgt)));
end

% prune out weak linkages
thrP = prctile(distClickE,p.pruneThr);

[values, rows, cols, isolated] = prune(tempN, distClickE, rows1, cols1, thrP);


values = real(values);
nodesRemaining = length(values);
if nodesRemaining < p.minClust
    fprintf('Too few nodes remaining after pruning, skipping iteration \n')
    spectraMean = [];%mean(specClickTf_norm);
    return % if there aren't nodes/edges left, skip iteration
end

% Run Gephi clustering routine
file2Write = sprintf('click_output.gv');
fid = fopen(file2Write, 'W');
fprintf(fid, 'strict graph {\n');
for vidx=1:length(values)
    fprintf(fid, 'n%d -- n%d [weight=%f, style=invis];\n', ...
        rows(vidx),cols(vidx), values(vidx));
end
fprintf(fid, '}\n');
fclose(fid);

disp('clustering')
[status,result] = system(sprintf('"%s" -Xms512m -Xmx20000m -Dfile.encoding=Cp1252 -classpath %s;%s ClusterGephi.ClusterMain "%s" "%g" "%g"',...
    javaPathVar,classPathVar,toolkitPath,file2Write,p.modular,p.pgThresh/100));
disp('done clustering')
if status~=0
    error(result)
end


% Read in output table from csv created by Gephi
[dataTable,labels] = xlsread(filename);
% eigenVec = dataTable(4,:)>=.5;
% pgRank = dataTable(3,:)>=prctile(dataTable(3,:),75);
clusters = dataTable(3,:)';
clustersAll = dataTable(1,:)';
nodeNums = char(labels(1,:)');
nodeNumsAll = char(labels(1,:)');
% rankVals = dataTable(3,pgRank);
nodeNums = str2double(cellstr(nodeNums(:,2:end)));
nodeNumsAll =  str2double(cellstr(nodeNumsAll(:,2:end)));
clustBins = 0:max(clusters);
counts = histc(clusters,clustBins);
keepClust = find(counts >= p.minClust);
clustSizes = counts(keepClust);
clustNums = clustBins(keepClust);
spectraMean = [];
spectraHolder = {};
%spectraStd=[];
clickAssign = {};
clickAssignAll = {};

% Organize output into cell arrays by clusters
if ~isempty(keepClust)
    for i4 = 1:length(keepClust)
        linearSpec = 10.^(specClickTf_norm(nodeNums(clusters==clustNums(i4)),:)./20);
        spectraMean(i4,:) = 20*log10(mean(linearSpec));
        % spectraStd(i4,:) = std(specClickTf_norm(nodeNums(clusters==clustNums(i4)),:));
        clickAssignAll{i4} = nodeNumsAll(clustersAll == clustNums(i4));
        clickAssign{i4} = nodeNums(clusters==clustNums(i4));
        spectraHolder{i4} = specClickTf_norm(nodeNums(clusters==clustNums(i4)),:);
    end
end
% 
% if length(nodenums)>minClust
%     spectraMean = mean(specClickTf_norm(nodenums,:));
% 	spectraStd = std(specClickTf_norm(nodenums,:));
%     clickAssign = nodenums;
% 	spectraHolder = specClickTf_norm(nodenums,:);
% end

 1;