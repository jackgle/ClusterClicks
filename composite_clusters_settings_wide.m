function s = composite_clusters_settings_wide

% Settings for running composite_clusters.m
% Consider making separate versions to keep track of settings used on
% different projects

% File locations and names
s.siteName = 'JAX11D_PG0_PR98'; % used to name output files
s.inDir = 'D:\JAX11D\TPWS\ClusterBins_FPincluded'; 
s.inFile = 'JAX11D_disk08_Delphin_clusters_diff_PG0_PR98_MIN25_MOD0_PPmin120_FPincl';%'MC_GC_DT_01-02_autoCluster_90_2000ofEach.mat';
s.outDir = 'D:\JAX11D\TPWS\ClusterBins_demotest';

% Gephi and Java paths:
s.javaPathVar = 'C:\Program Files\Java\jre6\bin\java.exe';
s.classPathVar = 'E:\workspace\ClusterGephi_sio\bin';
s.toolkitPath = 'E:\workspace\ClusterGephi_sio\gephi-toolkit-0.8.7-all\gephi-toolkit.jar';

% Saving?
s.saveOutput = 1; %set to 1 to save output file and figs, else 0

%%%% Similarity %%%%
% choose if you want to include ICI **OR** click rate in similarity calculation
s.iciModeTF = 1; % 1 if you want to use modal ICI (time between clicks)
%OR
s.iciDistTF = 0;% 1 if you want to compare ici distributions
%OR
s.cRateTF = 0; % 1 if you want to use click rate (# of clicks per second)
%%%
s.correctForSaturation = 1; % 1 if you want to look for minor ICI peaks in 
% cases where clicking is so dense that individual ICIs are obscured. This
% helps with dolphins, but may hurt if you are trying to pull out ships and
% rain too. Only for modal ICI
s.specDiffTF = 1; % set to 1 to use spectral 1st derivatives for correlation

%%%% Distribution Pruning %%%%
s.stIdx = 5;
s.edIdx = 190;
s.maxICI = 61;
s.minICI = 1;

%%%% Clustering %%%%
s.minClust = 10; % minimum number of bins required for a cluster to be retained.
s.pruneThr = 95; % Percentage of edges between nodes that you want to prune.
s.pgThresh = 0; % Percentile of nodes to remove from network using PageRank weights.
% e.g. If you use 25, nodes with PR in the lowest 25th percentile will be
% pruned out.
s.modular = 0; % If you use a number other than 0, modularity algorithm will
% be used instead of chinese whispers. Not recommended.
s.maxClust = 10000;% maximum number of bins to cluster. If you have more than
% this, a random subset of this size will be selected.
s.subSampOnlyOnce = 1; % if your input contains more than maxClust clicks, they 
% will be subsampled. If subSampOnlyOnce = 1, then a subsample will
% be selected, and it will be reclustered N times. This ends up looking at
% fewer clicks, but the best set of final clusters isn't chosen based on
% the simplest subset. It's also faster.
% If subSampOnlyOnce = 0, then a new subsample will be selected on each of 
% N iterations. This looks at more signals, but risks that the final
% clusters will be chosen from the subset that happened to have the least
% variability.
s.minClicks = 50; % minimum number of clicks per bin that you want to consider
% higher number makes cleaner clusters, but may miss things that click
% slowly.

% Number of clusterings to use for evidence accumulation
s.N = 1; % bigger is theoretically more robust, but takes longer

%%%% Plotting %%%%
s.subPlotSet = 1; % Set to 1 if you want plots with each click type as a subplot
s.indivPlots = 0; % Set to 1 if you want separate plots for each click type

