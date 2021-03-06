% cluster_bins
% Cluster clicks in n-minute bins by spectral similarities using Gephi API
%%% INPUT:
% Directory of TPWS files. These need to include frequency vector
% associated with spectra, or add your freq vector to script.
%%% OUTPUT:
% One (large) mat file containing dominant click type(s) in each N minute bin
% during the period spanned by the TPWS files. Includes
%   sumSpec   - mean spectrum/spectra for each bin,
%   nSpec     - # of clicks in bin associated with each summary spec,
%   percSpec  - percentage of clicks in bin associated with each summary spec,
%   clickRate - click rate distribution(s) per bin,
%   dTT       - ICI distribution(s) per bin,
%   cInt      - # clicks per bin,
%   tInt      - bin start and end times,
% Also includes setting values described below.

%%% NOTE: Output is saved after each TPWS file is processed in case of crash,
% but data is stored cumulatively, so if all goes well, you only need the
% last (biggest) file.
% kef 10/14/2016
%%
clearvars

%%% Java/Gephi toolkit setup
javaPathVar = 'C:\Program Files\Java\jre6\bin\java.exe';
classPathVar = ' E:\workspace\ClusterGephi_sio\bin';
toolkitPath = 'E:\workspace\ClusterGephi_sio\gephi-toolkit-0.8.7-all\gephi-toolkit.jar';

%%% set inputs and setting values
siteName = 'GofMX_MP03'; % First few letters of TPWS file names

% directory where those TPWS files live
inDir = 'G:\MP\MP03_TPWS';
outDir = 'G:\MP\MP03_TPWS\cluster_bins_output';
clusteredFileName = 'G:\MP\Matched\MP01-09_clusters_ALL_cell.mat';
cS = load(clusteredFileName,'tInt');
% pruneItrVals = [70,80,90,95,97,99];
% for iPrune = 1:length(pruneItrVals)
%%% Clustering parameter choices
p.minClust = 100; % minimum number of clicks required for a cluster to be retained.
% Think about how fast your species click, group sizes, and how many clicks they make
% per N minutes...
p.pruneThr = 97; % Percentage of edges between nodes that you want to prune.
% Pruning speeds up clustering, but can result in isolation (therefore
% loss) of rare click types.
p.pgThresh = 0; % Percentile of nodes to remove from metwork using PageRank weights.
% e.g. If you use 25, nodes with PR in the lowest 25th percentile will be
% pruned out.
p.modular = 0; % if you use a number other than 0, modularity algorithm will be used
% instead of chinese whispers for community detection. Not recommended.
% In the modularity algorithm, this parameter influences the number of
% communities detected. 1 = no bias, >1 bias toward fewer communities, <1,
% bias toward more communities.

p.falseRM = 1; % Want to remove false positives? Working on removing the
% need for manual false positive ID step.

%%% Frequencies you want to compare clicks across:
% comparing across the full bandwidth tends to reduce differences between click
% types. Results are typically better if you focus comparing the region
% where frequencies differ most.
p.stIdx = 1; % index of start freq in vector "f" from TPWS
p.edIdx = 121; % index of end freq

%%% Vectors to use for binning ICI and click rate
p.barInt = 0:.01:.5; % ICI bins in seconds (minICI:resolution:maxICI)
p.barRate = 1:1:30; % click rate in clicks per second (minRate:resolution:maxRate)
% Option to enforce a minimum recieved level (dB peak to peak), and only
% cluster high-amplitude clicks, which tend to have cleaner spectra.
p.ppThresh = 120;

%%% Time bins: clicks are clustered by time bin. How long of a bin do you
% want to consider?
% Larger bins -> more clicks, better representation for slow clicking species.
% But large time bins mean more click counts -> longer computation times,
% or subsampling
p.timeStep = 5; % bin duration in mins


%% Begin calculations
cd(inDir);

% initialize data holders
sumSpec = {};
nSpec = {};
percSpec = {};
cInt = [];
tInt = [];
dTT = {};
clickRate = {};
cIdx = 1;
nIsolated = {};
clusteredTF = [];
ttppNames = dir([siteName,'*_TPWS1.mat']);

fdNames = dir([siteName,'*_FD1.mat']);
fdAll = [];

% Load all FD files, in case for some reason they don't line up with TPWS
% files. This may be unnecessary now, and it's slower (input to setdiff)
falseStr = '_FPincl';

if p.falseRM
    for i0 = 1:length(fdNames)
        zFD = [];
        load(fdNames(i0).name,'zFD');
        fdAll = [fdAll;zFD];
    end
    falseStr = '_FPremov';
end

fkeep = [];
noCluster = 0;
for itr = 1:length(ttppNames)
    thisFile = ttppNames(itr).name;
    MTT = [];
    MPP = [];
    MSP = [];
    zFD = [];
    % FDname = strrep(thisFile,'TPWS','FD');
    load(thisFile,'MPP','MTT','MSP','f')
    if ~isempty(f)
        fkeep = f;
    elseif isempty(f) && ~isempty(fkeep)
        f = fkeep;
    else
        disp('Error: Missing frequency vector in input file.')
        break
    end
    % load(FDname,'zFD')
    % make output file name that incorporates settings:
    outName = strrep(thisFile,'TPWS1',sprintf('toofewtocluster_MIN%d',p.minClust));
    
    % remove false positive clicks
    [~, keepIdx] = setdiff(MTT,fdAll);
    MTT = MTT(keepIdx);
    MSP = MSP(keepIdx,:);
    MPP = MPP(keepIdx);
    
    % remove low amplitude clicks
    ppKeep = MPP>=p.ppThresh;
    MSP = MSP(ppKeep,:);
    MPP = MPP(ppKeep);
    MTT = MTT(ppKeep);
    
    if isempty(MTT) % go to next file if no clicks left in this one.
        continue
    end
    
    % Build vector of bin starts/ends.
    dateInterval = floor(MTT(1)):datenum([0,0,0,0,p.timeStep,0]):ceil(MTT(end));
    
    % Figure out which clicks are in which time bin.
    [testN,testBin] = histc(MTT,dateInterval);
    
    idxer = 1:length(testBin);
    % loop over bins
    for iC = 1:1:length(dateInterval)-1
        tTemp = [dateInterval(iC),dateInterval(iC+1)];
        nClicks = testN(iC); % store number of clicks in interval
        idSet = idxer(testBin == iC);
        
        if nClicks >0 && isempty(intersect(tTemp,cS.tInt,'rows'));
            specSet = MSP(idSet,:);
            ppSet = MPP(idSet);
            ttSet = MTT(idSet);
            
            nIsolated{cIdx,1} = NaN;
            sumSpec{cIdx,1} = calc_norm_spec_mean(specSet); % store summary spectra
            nSpec{cIdx,1} = nClicks; % store # of clicks associated with each summary spec
            percSpec{cIdx,1} = 1; % store % of clicks associated with each summary spec
            cInt(cIdx,1) = nClicks; % store number of clicks in interval
            tInt(cIdx,:) = [dateInterval(iC),dateInterval(iC+1)]; % store start and end time of interval
            dTT{cIdx,1} = histc(diff(sort(ttSet))*24*60*60,p.barInt)';
            clickRate{cIdx,1} = histc(1./(diff(sort(ttSet))*24*60*60),p.barRate)';
            clusteredTF(cIdx,:) = 0;
            cIdx = cIdx +1;
            noCluster = 0;
            
        end
        if mod(iC,200) == 0
            fprintf('done with bin # %d of %d, file %d\n',iC-1,length(dateInterval),itr)
        end
    end
    % save output
    save(fullfile(outDir,outName),'dTT','sumSpec','percSpec',...
        'cInt','tInt','p','f','nSpec', 'clickRate','nIsolated','clusteredTF')
    
    
end
