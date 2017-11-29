function [siteName, inDir, outDir, p] = cluster_bins_settings

%%% Set inputs and setting values for cluster_bins %%%

% Deployment info
siteName = 'MC02'; % First few letters of TPWS file names

% Folder info:
inDir = 'H:\MC\MC02_TPWS'; % where the TPWS files live, 
outDir = 'H:\MC\MC02_TPWS\demo'; % where to save outputs

% Java/Gephi toolkit setup
p.javaPathVar = 'C:\Program Files\Java\jre6\bin\java.exe';
p.classPathVar = ' E:\workspace\ClusterGephi_sio\bin';
p.toolkitPath = 'E:\workspace\ClusterGephi_sio\gephi-toolkit-0.8.7-all\gephi-toolkit.jar';

%%% Clustering parameter choices
p.minClust = 100;% minimum number of clicks required for a cluster to be retained.
% Think about how fast your species click, group sizes, and how many clicks they make
% per N minutes...
p.pruneThr = 95; % Percentage of edges between nodes that you want to prune.
% Pruning speeds up clustering, but can result in isolation (therefore
% loss) of rare click types.
p.pgThresh = 0; % Percentile of nodes to remove from network using PageRank weights.
% e.g. If you use 25, nodes with PR in the lowest 25th percentile will be
% pruned out.
p.modular = 0; % if you use a number other than 0, modularity algorithm will be used
% instead of chinese whispers for community detection. Not recommended.
% In the modularity algorithm, this parameter influences the number of
% communities detected. 1 = no bias, >1 bias toward fewer communities, <1,
% bias toward more communities.

p.plotFlag = 1; % Want plots? Turn this off for large jobs, but useful for
% seeing how your clusters are coming out when choosing parameters above.
p.falseRM = 0; % Want to remove false positives? Working on removing the 
% need for manual false positive ID step.

%%% Frequencies you want t6o compare clicks across:
% comparing across the full bandwidth tends to reduce differences between click
% types. Results are typically better if you focus comparing the region
% where frequencies differ most.
p.stIdx = 1; % index of start freq in vector "f" from TPWS
p.edIdx = 121; % index of end freq

%%% Vectors to use for binning ICI and click rate
p.barInt = 0:.01:.6; % ICI bins in seconds (minICI:resolution:maxICI)
p.barRate = 1:1:30; % click rate in clicks per second (minRate:resolution:maxRate)

p.diff = 1;% compare first derivative of spectra if 1. If 0 spectra will be compared normally.

% Option to enforce a minimum recieved level (dB peak to peak), and only
% cluster high-amplitude clicks, which tend to have cleaner spectra.
p.ppThresh = 120;

%%% Time bins: clicks are clustered by time bin. How long of a bin do you
% want to consider?
% Larger bins -> more clicks, better representation for slow clicking species.
% But large time bins mean more click counts -> longer computation times,
% or subsampling
p.timeStep = 5; % bin duration in mins
p.maxNetworkSz = 5000; % maximum # of clicks allowed in a network.
% If there are more clicks in a time bin than this number, a random subset
% will be selected for clustering. Your computer will need to handle
% maxNetworkSz^2 edges, so more RAM -> larger networks.
