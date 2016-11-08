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
siteName = 'GC'; % First few letters of TPWS file names

% directory where those TPWS files live
inDir = 'F:\GOM_clickTypePaper_detections\TPWS\GC01_02_03_TPWS';
outDir = 'F:\GOM_clickTypePaper_detections\TPWS\GC01_02_03_TPWS\Cluster2016';

%%% Clustering parameter choices
p.minClust = 100; % minimum number of clicks required for a cluster to be retained.
% Think about how fast your species click, group sizes, and how many clicks they make
% per N minutes...
p.pruneThr = 90; % Percentage of edges between nodes that you want to prune.
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

p.plotFlag = 1; % Want plots? Turn this off for large jobs, but useful for
% seeing how your clusters are coming out when choosing parameters above.

%%% Frequencies you want to compare clicks across:
% comparing across the full bandwidth tends to reduce differences between click
% types. Results are typically better if you focus comparing the region
% where frequencies differ most.
p.stIdx = 1; % index of start freq in vector "f" from TPWS
p.edIdx = 81; % index of end freq

%%% Vectors to use for binning ICI and click rate
p.barInt = 0:.01:.5; % ICI bins in seconds (minICI:resolution:maxICI)
p.barRate = 1:1:30; % click rate in clicks per second (minRate:resolution:maxRate)

p.diff = 0;% compare first derivative of spectra if 1. If 0 spectra will be compared normally.

% Option to enforce a minimum recieved level (dB peak to peak), and only
% cluster high-amplitude clicks, which tend to have cleaner spectra.
p.ppThresh = 110;

%%% Time bins: clicks are clustered by time bin. How long of a bin do you
% want to consider?
% Larger bins -> more clicks, better representation for slow clicking species.
% But large time bins mean more click counts -> longer computation times,
% or subsampling
p.timeStep = 5; % bin duration in mins
p.maxNetworkSz = 4000; % maximum # of clicks allowed in a network.
% If there are more clicks in a time bin than this number, a random subset
% will be selected for clustering. Your computer will need to handle
% maxNetworkSz^2 edges, so more RAM -> larger networks.



%% Begin calculations
cd(inDir);
if ~isdir(outDir) % check if outdir exists, make it if not.
    mkdir(outDir)
end

% initialize data holders
sumSpec = {};
nSpec = {};
percSpec = {};
cInt = [];
tInt = [];
dTT = {};
cIdx = 1;

ttppNames = dir([siteName,'*_TPWS1.mat']);

fdNames = dir([siteName,'*_FD1.mat']);
fdAll = [];

% Load all FD files, in case for some reason they don't line up with TPWS
% files. This may be unnecessary now, and it's slower (input to setdiff)
for i0 = 1:length(fdNames)
    zFD = [];
    load(fdNames(i0).name,'zFD');
    fdAll = [fdAll;zFD];
end
fkeep = [];
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
    if p.diff
        outName = strrep(thisFile,'TPWS1',sprintf('clusters_diff_PG%d_PR%d_MIN%d_MOD%d_PPmin%d_noFP',...
            p.pgThresh, p.pruneThr, p.minClust, p.modular,p.ppThresh));
    else
        outName = strrep(thisFile,'TPWS1',sprintf('clusters_PG%d_PR%d_MIN%d_MOD%d_noFP',...
            p.pgThresh, p.pruneThr, p.minClust, p.modular));
    end
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
        nClicks = testN(iC); % store number of clicks in interval

        if nClicks >= p.minClust*(1+(p.pgThresh/100))
            idSet = idxer(testBin == iC);
            specSet = MSP(idSet,:);
            ppSet = MPP(idSet);
            ttSet = MTT(idSet);
                        
            % if the interval contains more than 4000 good clicks, randomly
            % select a subset for clustering
            if nClicks > p.maxNetworkSz
                fprintf('%d clicks in bin, selecting subset of %d\n',nClicks,p.maxNetworkSz)
                rList = randperm(nClicks,p.maxNetworkSz);
                specSet = specSet(rList,:);
                ppSet = ppSet(rList);
                ttSet = ttSet(rList);
            end
            
            % Cluster
            [spectraMean,clickAssign,~,specHolder] = cluster_clicks(specSet,...
                p,javaPathVar,classPathVar,toolkitPath);
            
            
            % If you finish clustering with populated cluster(s)
            if ~isempty(clickAssign)
                % Calculate mean spectra, click rates, and
                % interclick-intervals
                % [~,idIdx,~] = intersect(zID(:,1),ttSet);
                % manLabels = zID(idIdx,2);
                sizeCA = zeros(size(clickAssign));
                dtt = zeros(size(clickAssign,2),length(p.barInt));
                cRate = zeros(size(clickAssign,2),length(p.barRate));
                for iS = 1:size(clickAssign,2)
                    sizeCA(1,iS) = size(clickAssign{iS},1);
                    dtt(iS,:) = histc(diff(sort(ttSet(clickAssign{iS})))...
                        *24*60*60,p.barInt);
                    cRate(iS,:) = histc(1./(diff(sort(ttSet(clickAssign{iS})))...
                        *24*60*60),p.barRate);
                end
                
                
                sumSpec{cIdx,1} = spectraMean; % store summary spectra
                nSpec{cIdx,1} = sizeCA; % store # of clicks associated with each summary spec
                percSpec{cIdx,1} = sizeCA./size(specSet,1); % store % of clicks associated with each summary spec
                cInt(cIdx,1) = nClicks; % store number of clicks in interval
                tInt(cIdx,:) = [dateInterval(iC),dateInterval(iC+1)]; % store start and end time of interval
                dTT{cIdx,1} = dtt;
                clickRate{cIdx,1} = cRate;
                cIdx = cIdx +1;
                %%
                if p.plotFlag
                    % plot mean spectra
                    figure(1);clf
                    subplot(2,2,1)
                    hold on
                    hR = rectangle('Position',[f(p.stIdx),0,f(p.edIdx)-f(p.stIdx),...
                        1],'FaceColor',[.99,.92,.8],'EdgeColor',[.99,.92,.8]);
                    plot(f,spectraMean)
                    ylim([min(min(spectraMean)),max(max(spectraMean))])
                    hold off
                    box on
                    xlim([10,70])
                    ylabel('Normalized Amplitude','Fontsize',12)
                    xlabel('Frequency (kHz)','Fontsize',12)
                    
                    % Plot ICI distribution
                    subplot(2,2,3)
                    plot(p.barInt,dtt)
                    xlim([0,.4])
                    ylabel('Counts','Fontsize',12)
                    xlabel('ICI (sec)','Fontsize',12)
                    
                    % Plot click rate distribution
                    subplot(2,2,4)
                    plot(p.barRate,cRate)
                    xlim([0,30])                   
                    ylabel('Counts','Fontsize',12)
                    xlabel('Click Rate (clicks/sec)','Fontsize',12)
                    
                    % plot concatonnated spectra. Need to handle case where
                    % multiple click types are present.
                    subplot(2,2,2)
                    allSpectra = cell2mat(specHolder');
                    imagesc(1:size(allSpectra,1),f,allSpectra')
                    delimLocs = cumsum(sizeCA);
                    % draw delimiters btwn click spectra if there are
                    % multiple types
                    if size(delimLocs,2)>1                            
                        hold on
                        for iDlim = 1:size(delimLocs,2)-1

                            plot([delimLocs(iDlim),delimLocs(iDlim)],...
                                [min(f),max(f)],'k','LineWidth',3)
                        end
                        hold off
                    end
                    set(gca,'ydir','normal','Fontsize',12)
                    ylabel('Frequency (kHz)','Fontsize',12)
                    xlabel('Click Number','Fontsize',12)
                    % pause
                    1;
                end
            end
        end
        if mod(iC,200) == 0
            fprintf('done with bin # %d of %d, file %d\n',iC-1,length(dateInterval),itr)
        end    
    end
    
    % save output
    save(fullfile(outDir,outName),'dTT','sumSpec','percSpec',...
        'cInt','tInt','p','f','nSpec', 'clickRate')
    
end
