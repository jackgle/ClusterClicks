% composite_clusters

% Clusters a set of mean spectra and ici distributions to identify major
% click types in a dataset.
% ** Run this on the output from cluster_bins.m **

%% Setup %%

clearvars

% Modify composite_clusters_settings.m to import your site- or project-
% specific preferences.
% You can save different versions of cluster_bins_settings, so you don't
% have to overwrite old settings. (TODO: add diary).

s = composite_clusters_settings;

colormap(jet)
%%
% Check for output directory
if ~exist(s.outDir,'dir')
    mkdir(s.outDir)
end
cd(s.outDir)

inFileName = fullfile(s.inDir,s.inFile);
% load data
load(inFileName)
inFileStub = s.siteName;

%%%%%%%%%% Begin main functionality %%%%%%%%%%%%
%% Normalize everything
% Spectra:
% Put click spectra into a matrix
if iscell(sumSpec)
    sumSpecMat = vertcat(sumSpec{:});
    nSpecMat = horzcat(nSpec{:})';
else
    sumSpecMat = sumSpec;
    nSpecMat = nSpec;
end
useBins = nSpecMat >= s.minClicks;
[specNorm,diffNormSpec] = spec_norm_diff(sumSpecMat(useBins,:),s.stIdx,s.edIdx);

% ICIs:
% move ICI distributions into a matrix and normalize.
if iscell(dTT)
    for iD = 1:size(dTT,1)
        % fixes some problems with past versions of inputs, 
        % eventually shouldn't need this, but it's here just to check that
        % things are oriented right for cell array to matrix conversion.
        if isempty(dTT{iD})
            dTT{iD} = zeros([1,length(p.barInt)]);
        end
        if size(dTT{iD},1)>size(dTT{iD},2)
            dTT{iD} = dTT{iD}';
        end
    end
    dTTmat = vertcat(dTT{:});
else
    dTTmat = dTT;
end
dTTmat = dTTmat(useBins,1:s.maxICI); % truncate if needed to ignore high ici peaks
dTTmatNorm1 = dTTmat./repmat(sum(dTTmat,2),1,size(dTTmat,2));
dTTmatNorm = dTTmatNorm1./repmat(max(dTTmatNorm1,[],2),1,size(dTTmat,2));

% Handle ICI distributions that have saturation issues due to overlapping
% clicking.
% find dTT rows where first bin is largest
[~,iciModeIdx] = max(dTTmatNorm(:,s.minICI:end),[],2);

% find secondary ICI peak in saturated ICI distributions
if s.correctForSaturation
    saturatedSet = find(iciModeIdx <= 3);
    
    for iSat = 1:length(saturatedSet)
        thisDTTsmoothDiff = diff(dTTmatNorm(saturatedSet(iSat),1:s.maxICI));
        minIdx = find(thisDTTsmoothDiff(2:end)>0, 1,'first')+1;
        if isempty(minIdx)
            minIdx = 1;
        end
        trucatedICIDist = dTTmatNorm(saturatedSet(iSat),minIdx:s.maxICI);
        trucatedICIDistNorm = trucatedICIDist/max(trucatedICIDist);

        [mVal,mTemp] = max(trucatedICIDistNorm,[],2);
        
        % compute rough peak prominence metric
        if mTemp>1 && mTemp<length(trucatedICIDistNorm)
            peakProm = mVal-((trucatedICIDistNorm(mTemp-1)+...
                trucatedICIDistNorm(mTemp+1))/2);
        elseif mTemp == 1
            peakProm = mVal-(trucatedICIDistNorm(mTemp+1)/2);
        elseif mTemp == length(trucatedICIDistNorm)
            peakProm = mVal-(trucatedICIDistNorm(mTemp-1)/2);
        end
        if peakProm>.2 && sum(trucatedICIDist)>.5% don't adjust unless the peak is strong enough
            iciModeIdx(saturatedSet(iSat)) = minIdx+mTemp-1;
        end
%         figure(1);clf;
%         plot(iciModeIdx(saturatedSet(iSat)),mVal,'*');
%         hold on;plot(dTTmatNorm(saturatedSet(iSat),:));plot([0,40],[0,0],'r');
%         plot(thisDTTsmoothDiff,'g');hold off
    end
end
iciModes = p.barInt(iciModeIdx+s.minICI) + p.barInt(2)./2;
% [iciDist,~,~] = ici_dist(dTTmatNorm);

% Click rates - convert to matrix and normalize
if iscell(clickRate)
    for iD = 1:size(clickRate,1)
        if isempty(clickRate{iD})
            clickRate{iD} = zeros([1,length(p.barRate)]);
        end
        if size(clickRate{iD},1)>size(clickRate{iD},2)
            clickRate{iD} = clickRate{iD}';
        end
    end
    cRateMat = vertcat(clickRate{:});
else
    cRateMat = clickRate;
end
cRateMat = cRateMat(useBins,:);
cRateNorm1 = cRateMat./repmat(sum(cRateMat,2),1,...
    size(cRateMat,2));
cRateNorm = cRateNorm1./repmat(max(cRateNorm1,[],2),1,size(cRateMat,2));
[~,cRateModeIdx] = max(cRateNorm,[],2);
cRateModes = p.barRate(cRateModeIdx) + p.barRate(2)./2;

%% Cluster N times for evidence accumulation/robustness
tempN = ceil(sqrt(size(cRateMat,1)*2));
% CoMat = zeros(tempN,tempN);
subSamp = 0; % flag automatically goes to true if subsampling.
isolatedSet = [];

for iEA = 1:s.N
    % Select random subset if needed
    if size(cRateMat,1)> s.maxClust
        subSamp = 1;
        excludedIn = sort(randperm(length(dTTmatNorm), s.maxClust));
        fprintf('Max number of bins exceeded. Selecting random subset of %d\n',s.maxClust)
        if s.subSampOnlyOnce
            % set flag back to zero if you only want to subsample once,
            % rather than taking a new subsample every time.
            subSamp = 0;
        end
    else
        excludedIn = 1:length(dTTmatNorm);
        subSamp = 0;
    end
    
    if subSamp || iEA == 1
        % Only do this on first iteration, or on every iteration if you are subsampling
        % find pairwise distances between spectra
        if s.specDiffTF
            [specDist,~,~] = spectra_dist(diffNormSpec(excludedIn,s.stIdx:s.edIdx));
        else
            [specDist,~,~] = spectra_dist(specNorm(excludedIn,s.stIdx:s.edIdx));
        end
        
        if s.iciModeTF % if true, use ici distributions for similarity calculations
            [iciModeDist,~,~,~] = ici_dist_mode(iciModes(excludedIn),p.barInt(s.maxICI));
            compDist = squareform(specDist.*iciModeDist,'tomatrix');
            disp('Clustering on modal ICI and spectral correlations')
        elseif s.iciDistTF
             % if true, use ici distributions for similarity calculations
            [iciDist,~,~] = ici_dist(dTTmatNorm(excludedIn,s.minICI:s.maxICI));
            compDist = squareform(specDist.*iciDist,'tomatrix');
            disp('Clustering on ICI distribution and spectral correlations')
        elseif s.cRateTF
            % use click rate distributions for similarity calculations
            [cRateDist,~,~] = ici_dist(cRateNorm(excludedIn,:));
            compDist = squareform(specDist.*cRateDist,'tomatrix');
            disp('Clustering on modal click rate and spectral correlations')
        else
            % if nothing, only use spectra
            compDist = squareform(specDist,'tomatrix');
            disp('Clustering on spectral correlations')
        end
        [gv_file,isolated] = write_gephi_input(compDist, s.minClust,s.pruneThr);
        
    end
    inputSet{iEA} = setdiff(excludedIn,isolated);
    fprintf('Clustering for evidence accumulation: %d of %d\n',iEA,s.N)
    
    % cluster
    nodeN = size(compDist,1);
    [nodeAssign,excludedOut,rankAssign] = run_gephi_clusters(nodeN,...
        s.minClust,s.modular,s.pgThresh,s.javaPathVar,s.classPathVar,s.toolkitPath,gv_file);
    
    % Recover from random permutation
    iL = 1;
    nodeSet = {};
    for iA = 1:length(nodeAssign)
        nodeSet{iA} = excludedIn(nodeAssign{iL});
        % update co-association matrix (Fred and Jain 2005)
        % CoMat(nodeSet{iA},nodeSet{iA}) = CoMat(nodeSet{iA},nodeSet{iA})+ 1/N;
        iL = iL+1;
    end
    nList{iEA} = excludedIn;
    ka(iEA,1) = length(nodeAssign);
    naiItr{iEA} = nodeSet;
    fprintf('found %d clusters\n',length(nodeAssign))
    isolatedSet(iEA,1) = length(setdiff(excludedIn,horzcat(nodeSet{:})));
end

% Best of K partitions based on NMI filkov and Skiena 2004
% Compute NMI
fprintf('Calculating NMI\n') 
[NMIList] = compute_NMI(nList,ka,naiItr,inputSet);
% the one with the highest mean score is the best
[bokVal,bokIdx] = max(sum(NMIList)./(size(NMIList,1)-1)); % account for empty diagonal.

nodeSet = naiItr{bokIdx};
% % Final cluster using Co-assoc mat.
% disp('Clustering using accumulated evidence')
% compDist2 = compDist;
% compDist2(find(CoMat<1))=0;
% [nodeSet,excludedOut,rankAssign] = cluster_nodes(compDist2,...
%        minClust,pruneThr,modular,pgThresh);

%% calculate means, percentiles, std devs
for iF = 1:length(nodeSet)
    % compute mean of spectra in linear space
    linearSpec = 10.^(specNorm(nodeSet{iF},:)./20);
    spectraMeanSet(iF,:)=  20*log10(nanmean(linearSpec));
    specPrctile{iF} = prctile(specNorm(nodeSet{iF},:),[25,75]);
    iciMean(iF,:) = nanmean(dTTmatNorm(nodeSet{iF},:));
    iciStd(iF,:) = nanstd(dTTmatNorm(nodeSet{iF},:));
    cRateMean(iF,:) = nanmean(cRateNorm(nodeSet{iF},:));
    cRateStd(iF,:) = nanstd(cRateNorm(nodeSet{iF},:));
end


%% plotting
barAdj = .5*mode(diff(p.barInt));
if s.subPlotSet
    n1 = 3; % number of rows of subplots, one subplot per type
    m1 = ceil(length(nodeSet)/n1); % number of columns of subplots
    figure(40);clf;figure(90);clf
    figure(80);clf;figure(81);clf;
    
    typeCount = size(nodeSet,2);
    for iF = 1:length(nodeSet)
        figure(40) % plot spectra means and percentiles
        subplot(n1,m1,iF)
        hold on
        plot(f,spectraMeanSet(iF,:),'-k','lineWidth',2)
        xlim([min(f),max(f)])
        legend(num2str(size(nodeSet{iF},2)),'location','Southeast')
        plot(f,specPrctile{iF},'--k','lineWidth',2)
        grid on
        
        figure(90) % plot spectra as heatmap
        subplot(n1,m1,iF)
        imagesc(1:length(nodeSet{iF}),f,specNorm(nodeSet{iF},:)')
        set(gca,'ydir','normal')
        ylim([min(f),max(f)])
        
        figure(80) % plot ICI distributions
        subplot(n1,m1,iF)
        errorbar(p.barInt(1:s.maxICI)+barAdj,iciMean(iF,:),zeros(size(iciStd(iF,:))),iciStd(iF,:),'.k')
        hold on
        hbI = bar(p.barInt(1:s.maxICI)+barAdj,iciMean(iF,:),1);
        xlim([0,p.barInt(s.maxICI)])
        ylim([0,1])
        
        figure(81) % plot click rate distributions
        subplot(n1,m1,iF)
        errorbar(p.barRate,cRateMean(iF,:),zeros(size(cRateStd(iF,:))),cRateStd(iF,:),'.k')
        hold on
        hbR = bar(p.barRate,cRateMean(iF,:),1);
        xlim([0,max(p.barRate)])
        ylim([0,1])
    end
    if s.saveOutput
        figName80 = fullfile(s.outDir,sprintf('%s_autoTypes_allICI',inFileStub));
        set(80,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25  10  7.5])
        figure(80)
        x80 = mxlabel(gcf,'ICI (sec)','FontSize',16);
        y80 = mylabel(gcf,'Normalized Counts','FontSize',16);
        print(80,'-dtiff','-r600',[figName80,'.tiff'])
        saveas(80,[figName80,'.fig'])
        
        figName81 = fullfile(s.outDir,sprintf('%s_autoTypes_allcRate',inFileStub));
        set(81,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25  10  7.5])
        figure(81)
        x81 = mxlabel(gcf,'Click Rate (clicks/sec)','FontSize',16);
        y81 = mylabel(gcf,'Normalized Counts','FontSize',16);
        print(81,'-dtiff','-r600',[figName81,'.tiff'])
        saveas(81,[figName81,'.fig'])
        
        set(40,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25  10  7.5])
        figName40 = fullfile(s.outDir,sprintf('%s_autoTypes_allMeanSpec',inFileStub));
        figure(40)
        x40 = mxlabel(gcf,'Frequency (kHz)','FontSize',16);
        y40 = mylabel(gcf,'Normalized Amplitude','FontSize',16);
        print(40,'-dtiff','-r600',[figName40,'.tiff'])
        saveas(40,[figName40,'.fig'])
        
        set(90,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25  10  7.5])
        figName90 = fullfile(s.outDir,sprintf('%s_autoTypes_allCatSpec',inFileStub));
        figure(90)
        x90 = mxlabel(gcf,'Click Number (kHz)','FontSize',16);
        y90 = mylabel(gcf,'Frequency (kHz)','FontSize',16);
        print(90,'-dtiff','-r600',[figName90,'.tiff'])
        saveas(90,[figName90,'.fig'])
    end
end

%% Individual plots
if s.indivPlots
    siteNameNo_ = strrep(s.siteName,'_','\_');
    typeCount = size(nodeSet,2);
    % close all
    for iF = 1:length(nodeSet)
        figure(400);clf
        
        subplot(1,4,1)
        errorbar(p.barInt(1:s.maxICI)+barAdj,iciMean(iF,:),zeros(size(iciStd(iF,:))),iciStd(iF,:),'.k')
        xlabel('ICI (sec)','FontSize',12)
        hold on
        hbI2 = bar(p.barInt(1:s.maxICI)+barAdj,iciMean(iF,:),1);
        xlim([0,p.barInt(s.maxICI)])
        ylim([0,1])
        ylabel('Relative Counts','FontSize',12)
        annotation('textbox',[.02 .82 .1 .1],'units','normalized','String',...
            siteNameNo_,'LineStyle','none','FontSize',14,'FontWeight','demi')
        set(gca,'box','on','fontsize',11)
        hold off
        ylim([0,1])
        
        subplot(1,4,2)
        annotation('textbox',[.22 .82 .1 .1],'units','normalized','String',...
            ['n = ',num2str(size(nodeSet{iF},2))],'LineStyle','none','FontSize',11)
        errorbar(p.barRate,cRateMean(iF,:),zeros(size(cRateStd(iF,:))),cRateStd(iF,:),'.k')
        hold on
        hbR2 = bar(p.barRate,cRateMean(iF,:),1);
        xlabel('Click Rate (clicks/sec)','FontSize',12)
        set(gca,'box','on','fontsize',11)
        hold off
        ylim([0,1])
        xlim([min(p.barRate),max(p.barRate)])
        
        hs3 = subplot(1,4,3);
        hold on
        hp2 = plot(f,spectraMeanSet(iF,:),'-k','lineWidth',2);
        set(gca,'box','on','FontSize',11)
        ylim([0,1])
        xlim([min(f),max(f)])
        plot(f,specPrctile{iF},':k','lineWidth',2)
        ylabel('Relative Amplitude','FontSize',12)
        xlabel('Frequency (kHz)','FontSize',12)
        hold off
        grid on
        set(gca,'XMinorTick', 'on')
        hs3Pos = get(hs3,'Position');
        
        subplot(1,4,4)
        imagesc(1:length(nodeSet{iF}),f,specNorm(nodeSet{iF},:)') ;
        h3im = gca;
        xlabel('Spectrum Number','FontSize',12)
        ylabel('Frequency (kHz)','FontSize',12)
        oldPos = get(h3im,'Position');
        set(h3im,'FontSize', 11)
        % ylim([f(stIdx),f(edIdx)])
        hc = colorbar;
        set(h3im,'Position',[oldPos(1),0.1404,hs3Pos(3),0.7846])
        caxis([0 1])
        
        set(hc,'Position',get(hc,'Position')-[0.01,0,.01,0],'Ylim',[0,1],...
            'FontSize',12)
        hcy = ylabel(hc,'Relative Amplitude','Rotation',-90);
        set(hcy,'Units','normalized')
        set(hcy,'Position',get(hcy,'Position')+[1,0,0],'Units','normalized')
        set(gca,'ydir','normal')
        set(gcf,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25 11 4])
        if s.saveOutput
            figName2 = fullfile(s.outDir,sprintf('%s_AutoType%d_longICI',inFileStub,iF));
            print(gcf,'-dtiff','-r600',[figName2,'.tif'])
            
            saveas(gcf,[figName2,'.fig'])
        end
    end
end
%%

for iTF = 1:length(nodeSet)
    Tfinal{iTF,1} = sumSpecMat(nodeSet{iTF},:);% all mean spectra
    Tfinal{iTF,2} = dTTmatNorm(nodeSet{iTF},:);% ici ditributions
    Tfinal{iTF,3} = diffNormSpec(nodeSet{iTF},:); % 1st diff spectra
    Tfinal{iTF,4} = iciModes(nodeSet{iTF}); % modal ICI vals
    Tfinal{iTF,5} =  spectraMeanSet(iTF,:); % mean spectrum
end
if s.saveOutput
    save(fullfile(s.outDir,[inFileStub,'_typesHR']),'Tfinal','nodeSet','NMIList',...
        'p','barAdj','f','nList','ka','naiItr','isolated','s','cInt','dTT',...
        'sumSpec','tInt', 'percSpec','nSpec','clickRate','specNorm',...
        'dTTmatNorm','cRateNorm','inFileName')
end