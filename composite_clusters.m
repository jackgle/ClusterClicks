% composite_clusters

% Clusters a set of mean spectra and ici distributions to identify major
% click types in a dataset. 
% ** Run this on the output from cluster_bins **

clearvars

siteName = 'MC01-03';
inDir = 'F:\GOM_clickTypePaper_detections\TPWS\MC01_02_03_TPWS\Cluster2016\';
inFile = 'MC03_disk14_clusters_PG0_PR95_MIN100_MOD0_noFP.mat';
outDir = 'F:\GOM_clickTypePaper_detections\TPWS\ClusterOct2016\';

%%%% Similarity %%%%
% choose if you want to include ICI **OR** click rate in similarity caculation
iciTF = 1; % 1 if you want to use ICI (time between clicks) 
cRateTF = 0; % 1 if you want to use click rate (# of clicks per second)

%%%% Distribution Pruning %%%%
stIdx = 1;
edIdx = 81;
maxICI = 41;
minICI = 4;

%%%% Clustering %%%%
minClust = 25; % minimum number of bins required for a cluster to be retained.
pruneThr = 95; % Percentage of edges between nodes that you want to prune.
pgThresh = 0; % Percentile of nodes to remove from metwork using PageRank weights.
% e.g. If you use 25, nodes with PR in the lowest 25th percentile will be
% pruned out.
modular = 0; % If you use a number other than 0, modularity algorithm will 
% be used instead of chinese whispers. Not recommended.
maxClust = 12000;% maximum number of bins to cluster. If you have more than
% this, a random subset of this size will be selected.

% Number of clusterings to use for evidence accumulation
N = 10; % bigger is theoretically more robust, but takes longer

%%%% Plotting %%%%
subPlotSet = 1; % Set to 1 if you want plots with each click type as a subplot
indivPlots = 0; % Set to 1 if you want separate plots for each click type

%% 
% Check for output directory
if ~exist(outDir,'dir')
    mkdir(outDir)
end
cd(outDir)

inFileName = fullfile(inDir,inFile);
% load data
load(inFileName)

%% Begin main functionality

% Put click spectra into a matrix
sumSpecMat = vertcat(sumSpec{:});
% find pairwise distances between spectra
[specDist,~,~,normSpec,diffNormSpec] = spectra_dist(sumSpecMat,stIdx,edIdx);

%%% Normalize ICIs
% move ICI distributions into a matrix
dTTmat = vertcat(dTT{:});
dTTmat = dTTmat(:,1:maxICI); % truncate if needed to ignore high ici peaks
% normalize ICI distributions
dTTmatNorm1 = dTTmat./repmat(sum(dTTmat,2),1,size(dTTmat,2));
dTTmatNorm = dTTmatNorm1./repmat(max(dTTmatNorm1,[],2),1,size(dTTmat,2));
% [iciDist,~,~] = ici_dist(dTTmatNorm);

%%% Normalize click rates
cRateMat = cell2mat(clickRate);
cRateNorm1 = cRateMat./repmat(sum(cRateMat,2),1,size(cRateMat,2));
cRateNorm = cRateNorm1./repmat(max(cRateNorm1,[],2),1,size(cRateMat,2));

if iciTF % if true, use ici distributions for similarity calculations
    [iciDist,~,~,iciModes] = ici_dist_mode(dTTmatNorm, p.barInt,minICI);    
    compDist = squareform(specDist.*iciDist,'tomatrix');
    disp('Clustering on modal ICI and spectral correlations')
elseif cRateTF 
    % use click rate distributions for similarity calculations
    [cRateDist,~,~,iciModes] = ici_dist_mode(cRateNorm, p.barRate,1);
    compDist = squareform(specDist.*cRateDist,'tomatrix');
    disp('Clustering on modal click rate and spectral correlations')
else
    % if nothing, only use spectra
    compDist = squareform(specDist,'tomatrix');
    disp('Clustering on spectral correlations') 
end

%% Cluster N times for evidence accumulation/robustness
tempN = ceil(sqrt(length(specDist)*2));
CoMat = zeros(tempN,tempN);
for iEA = 1:N
    fprintf('Clustering for evidence accumulation: %d of %d\n',iEA,N)
   
    % Select random subset if needed
    if size(compDist,1)>maxClust
        excludedIn = sort(randperm(length(dTTmatNorm),maxClust));
        fprintf('Max number of bins exceeded. Selecting random subset of %d\n',maxClust)
    else
        excludedIn = 1:length(dTTmatNorm);
    end
    
    % cluster
    [nodeAssign,excludedOut,rankAssign] = cluster_nodes(compDist(excludedIn,excludedIn),...
        minClust,pruneThr,modular,pgThresh);
    
    % Recover from random permutation
    iL = 1;
    nodeSet = {};
    for iA = 1:length(nodeAssign)
        nodeSet{iA} = excludedIn(nodeAssign{iL});
        % update co-association matrix (Fred and Jain 2005)
        CoMat(nodeSet{iA},nodeSet{iA}) = CoMat(nodeSet{iA},nodeSet{iA})+ 1/N;
        iL = iL+1;
    end
    nList{iEA} = excludedIn;
    ka(iEA,1) = length(nodeAssign);
    naiItr{iEA} = nodeSet;
    
   
end

% Best of K partitions based on NMI filkov and Skiena 2004
% Compute NMI
[NMIList] = compute_NMI(nList,ka,naiItr);
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
    linearSpec = 10.^(normSpec(nodeSet{iF},:)./20);
    spectraMeanSet(iF,:)=  20*log10(nanmean(linearSpec));
    specPrctile{iF} = prctile(normSpec(nodeSet{iF},:),[25,75]);
    iciMean(iF,:) = nanmean(dTTmatNorm(nodeSet{iF},:));
    iciStd(iF,:) = nanstd(dTTmatNorm(nodeSet{iF},:));
    cRateMean(iF,:) = nanmean(cRateNorm(nodeSet{iF},:));
    cRateStd(iF,:) = nanstd(cRateNorm(nodeSet{iF},:));
end


%% plotting
barAdj = .5*mode(diff(p.barInt));
if subPlotSet
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
        xlim([min(f),65])
        legend(num2str(size(nodeSet{iF},2)),'location','Southeast')
        plot(f,specPrctile{iF},'--k','lineWidth',2)
        grid on
        
        figure(90) % plot spectra as heatmap
        subplot(n1,m1,iF)
        imagesc(1:length(nodeSet{iF}),f,normSpec(nodeSet{iF},:)')
        set(gca,'ydir','normal')
        ylim([min(f),65])
        
        figure(80) % plot ICI distributions
        subplot(n1,m1,iF)
        errorbar(p.barInt(1:maxICI)+barAdj,iciMean(iF,:),zeros(size(iciStd(iF,:))),iciStd(iF,:),'.k')
        hold on
        hbI = bar(p.barInt(1:maxICI)+barAdj,iciMean(iF,:),1);
        xlim([0,p.barInt(maxICI)])
        ylim([0,1])

        figure(81) % plot click rate distributions
        subplot(n1,m1,iF)
        errorbar(p.barRate,cRateMean(iF,:),zeros(size(cRateStd(iF,:))),cRateStd(iF,:),'.k')
        hold on
        hbR = bar(p.barRate,cRateMean(iF,:),1);
        xlim([0,max(p.barRate)])
        ylim([0,1])
end

    figName80 = fullfile(outDir,sprintf('%s_autoTypes_allICI',siteName));
    set(80,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25  10  7.5])
    figure(80)
    x80 = mxlabel(gcf,'ICI (sec)','FontSize',16);
    y80 = mylabel(gcf,'Normalized Counts','FontSize',16);
    print(80,'-dtiff','-r600',[figName80,'.tiff'])
    saveas(80,[figName80,'.fig'])
    
    figName81 = fullfile(outDir,sprintf('%s_autoTypes_allcRate',siteName));
    set(81,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25  10  7.5])
    figure(81)
    x81 = mxlabel(gcf,'Click Rate (clicks/sec)','FontSize',16);
    y81 = mylabel(gcf,'Normalized Counts','FontSize',16);
    print(81,'-dtiff','-r600',[figName80,'.tiff'])
    saveas(81,[figName80,'.fig'])

    set(40,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25  10  7.5])
    figName40 = fullfile(outDir,sprintf('%s_autoTypes_allMeanSpec',siteName));
    figure(40)
    x40 = mxlabel(gcf,'Frequency (kHz)','FontSize',16);
    y40 = mylabel(gcf,'Normalized Amplitude','FontSize',16);
    print(40,'-dtiff','-r600',[figName40,'.tiff'])
    saveas(40,[figName40,'.fig'])
    
    set(90,'units','inches','PaperPositionMode','auto','OuterPosition',[0.25 0.25  10  7.5])
    figName90 = fullfile(outDir,sprintf('%s_autoTypes_allCatSpec',siteName));
    figure(90)
    x90 = mxlabel(gcf,'Click Number (kHz)','FontSize',16);
    y90 = mylabel(gcf,'Frequency (kHz)','FontSize',16);
    print(90,'-dtiff','-r600',[figName90,'.tiff'])
    saveas(90,[figName90,'.fig'])
end

%% Individual plots
if indivPlots
    siteNameNo_ = strrep(siteName,'_','\_');
    typeCount = size(nodeSet,2);
    % close all
    for iF = 1:length(nodeSet)
        figure(400);clf        
        
        subplot(1,4,1)
        errorbar(p.barInt(1:maxICI)+barAdj,iciMean(iF,:),zeros(size(iciStd(iF,:))),iciStd(iF,:),'.k')
        xlabel('ICI (sec)','FontSize',12)
        hold on
        hbI2 = bar(p.barInt(1:maxICI)+barAdj,iciMean(iF,:),1);
        xlim([0,p.barInt(maxICI)])
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
        imagesc(1:length(nodeSet{iF}),f,normSpec(nodeSet{iF},:)') ;
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
        figName2 = fullfile(outDir,sprintf('%s_AutoType%d_longICI',siteName,iF));
        print(gcf,'-dtiff','-r600',[figName2,'.tif'])
        saveas(gcf,[figName2,'.fig'])
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
[~,inFileStub,~] = fileparts(inFileName);
save(fullfile(outDir,[inFileStub,'_typesHR']),'Tfinal','nodeSet','NMIList',...
'p','maxClust','minClust','maxICI','barAdj','f','siteName','nList','ka','naiItr',...
'excludedOut','CoMat','iciTF','cRateTF','stIdx','edIdx','maxICI','minICI',...
'pruneThr','pgThresh','modular','N','cInt','dTT','sumSpec','tInt',...
'percSpec','nSpec','clickRate','normSpec','dTTmatNorm','cRateNorm','inFileName')
