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
%% Setup

clearvars

% modify cluster_bins_settings.m to import your site-specific preferences.
% you can save different versions of cluster_bins_settings, so you don't
% have to overwrite old settings. (TODO: add diary).
[siteName, inDir, outDir, p] = cluster_bins_settings;

colormap(jet)

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
    if p.diff
        outName = strrep(thisFile,'TPWS1',sprintf('clusters_diff_PG%d_PR%d_MIN%d_MOD%d_PPmin%d%s',...
            p.pgThresh, p.pruneThr, p.minClust, p.modular,p.ppThresh,falseStr));
    else
        outName = strrep(thisFile,'TPWS1',sprintf('clusters_PG%d_PR%d_MIN%d_MOD%d%s_PPmin%d%s',...
            p.pgThresh, p.pruneThr, p.minClust, p.modular,p.ppThresh,falseStr));
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
        idSet = idxer(testBin == iC);
        specSet = MSP(idSet,:);
        ppSet = MPP(idSet);
        ttSet = MTT(idSet);
        
        if nClicks >= p.minClust*(1+(p.pgThresh/100))

                        
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
            [spectraMean,clickAssign,~,specHolder,isolated] = cluster_clicks(specSet,...
                p,p.javaPathVar,p.classPathVar,p.toolkitPath);
            
            % If you finish clustering with populated cluster(s)
            if ~isempty(clickAssign)
                % Calculate mean spectra, click rates, and
                % interclick-intervals
                % [~,idIdx,~] = intersect(zID(:,1),ttSet);
                % manLabels = zID(idIdx,2);
                sizeCA = zeros(size(clickAssign));
                dtt = zeros(size(clickAssign,2),length(p.barInt));
                cRate = zeros(size(clickAssign,2),length(p.barRate));
                isolatedSet = zeros(size(clickAssign));
                meanSim = zeros(size(clickAssign));
                for iS = 1:size(clickAssign,2)
                    sizeCA(1,iS) = size(clickAssign{iS},1);
                    dtt(iS,:) = histc(diff(sort(ttSet(clickAssign{iS})))...
                        *24*60*60,p.barInt)';
                    cRate(iS,:) = histc(1./(diff(sort(ttSet(clickAssign{iS})))...
                        *24*60*60),p.barRate)';
                    isolatedSet(1,iS) = length(isolated);
                    % temporay calculation of within cluster similarity.
                    % disable soon!
%                     distClickTemp = pdist(specHolder{iS},'correlation');
%                     meanSim(1,iS) = mean(exp(-distClickTemp));
                end
                
                nIsolated{cIdx,1} = isolatedSet;
                sumSpec{cIdx,1} = spectraMean; % store summary spectra
                nSpec{cIdx,1} = sizeCA; % store # of clicks associated with each summary spec
                percSpec{cIdx,1} = sizeCA./size(specSet,1); % store % of clicks associated with each summary spec
                cInt(cIdx,1) = nClicks; % store number of clicks in interval
                tInt(cIdx,:) = [dateInterval(iC),dateInterval(iC+1)]; % store start and end time of interval
                dTT{cIdx,1} = dtt;
                clickRate{cIdx,1} = cRate;
                clusteredTF(cIdx,:) = 1;
                cIdx = cIdx +1;
                %meanSimilarity{cIdx,:} = meanSim; % disable soon!
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
                    xlim([f(p.stIdx),f(p.edIdx)])
                    ylabel('Normalized Amplitude','Fontsize',12)
                    xlabel('Frequency (kHz)','Fontsize',12)
                    
                    % Plot ICI distribution
                    subplot(2,2,3)
                    plot(p.barInt,dtt)
                    xlim([min(p.barInt),max(p.barInt)])
                    ylabel('Counts','Fontsize',12)
                    xlabel('ICI (sec)','Fontsize',12)
                    
                    % Plot click rate distribution
                    subplot(2,2,4)
                    plot(p.barRate,cRate)
                    xlim([min(p.barRate),max(p.barRate)])                   
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
                    drawnow
                    1;
                end
            else
                noCluster = 1;
            end
        elseif nClicks>0
            noCluster = 1;
        end
        if noCluster
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

%end