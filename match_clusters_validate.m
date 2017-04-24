clearvars
siteStr = 'MP3';
savDir = 'E:\Data\Papers\ClickClass2015\revamp';
fileToMatch = 'G:\MP\MP04_TPWS\cluster_bins_output\GofMX_MP04_disk09_Delphin_clusters_diff_PG0_PR97_MIN100_MOD0_PPmin120_FPremov.mat';
typeFile = 'F:\GOM_clickTypePaper_detections\TPWS\MC_GC_DTmerge\iciFix\MC_GC_DT_MP_DC_autoCluster_95_1000ofEach_1-121_diff_min200_icifix_typesHR.mat';
% siteStr = 'HAT_05_A';
% savDir = 'E:\Data\John Reports\HAT05A\dolphins';
% fileToMatch = 'G:\MC\Matched\MC01-11_clusters_ALL_cell.mat';
% typeFile = 'E:\Data\John Reports\HAT05A\dolphins\HAT_A_05_subsetof8000_matrix_typesHR.mat';

manualValidate = 1; % set to 1 if you want to manually assign an ID to bins with no good matches.
% else make it 0

load(fileToMatch);
load(typeFile,'Tfinal','iLine');
Tfinal = Tfinal(:,:);

% Prune Tfinal if needed
p.edIdx = 121;
nTemplates = size(Tfinal,1);
cMax = 0.25;
visualize = 0;
minICI = 1;
maxICI = 41;
% Check for output directory

% check 100 entries
myInt = 56; %unique(round(1:size(tInt,1)/200:size(tInt,1)));
if ~exist(savDir,'dir')
    mkdir(savDir)
end

autoScore = cell(size(tInt,1),1);
autoID = cell(size(tInt,1),nTemplates);
labeledCountsAll = zeros(size(tInt,1),nTemplates+1);
pSpecAll = cell(size(tInt,1),1);
cAll = zeros(size(tInt,1),1);
manVAuto = [];
%%
specCentroids = [];
iciCentroids = [];%nan(size(Tfinal,1),1);
nTf = ones(size(Tfinal,1),1);
normTfinalSpec = [];

for iTF = 1:size(Tfinal,1)
    iciCentroids(iTF,:) = mean(Tfinal{iTF,2});
    specCentroid1 = mean(Tfinal{iTF,1}(:,p.stIdx:p.edIdx));
    specCentroids(iTF,:) = (specCentroid1-min(specCentroid1))./max(specCentroid1-min(specCentroid1));
    nTf(iTF) = size(Tfinal{iTF},1);
    normTfinalSpec1 = Tfinal{iTF,1}(:,p.stIdx:p.edIdx);
    normTfinalSpec1 = normTfinalSpec1-repmat(min(normTfinalSpec1,[],2),1,length(p.stIdx:p.edIdx));
    normTfinalSpec{iTF,1} = normTfinalSpec1./repmat(max(normTfinalSpec1,[],2),1,length(p.stIdx:p.edIdx));
end

minBinIdx = 1;
if manualValidate
    figure(222);clf;
    subplot(1,2,1)
    plot(f(p.stIdx:p.edIdx),specCentroids)
    legend(num2str([1:nTemplates]'))
    xlim([10,65])
    subplot(1,2,2)
    plot(p.barInt(1:maxICI),iciCentroids)
end

for i0 = 1:length(sumSpec)
    % normalize mean spectra for this bin and compute diff
    if cInt(i0)>1
        minSSsection = min(sumSpec{i0},[],2);
        specClickTf_minNorm = (sumSpec{i0} - ...
            minSSsection(:,ones(1,size(sumSpec{i0},2))));
        maxSSsection = max(specClickTf_minNorm,[],2);
        specClickTf_norm = specClickTf_minNorm./maxSSsection(:,ones(1,size(specClickTf_minNorm,2)));
        specClickTf_norm_short = specClickTf_norm(:,p.stIdx:p.edIdx);
        %specClickTf_norm_short = diff(specClickTf_norm_short,1,2);
        
        % determine modal ici(s)
        if size(dTT{i0},1)>size(dTT{i0},2)
            dTT{i0} = dTT{i0}';
        end
        dTTmat = vertcat(dTT{i0});
        % normalize ICI distributions
        dTTmatNorm1 = dTTmat./repmat(sum(dTTmat,2),1,size(dTTmat,2));
        dTTmatNorm = dTTmatNorm1./repmat(max(dTTmatNorm1,[],2),1,size(dTTmat,2));
        
        [~,iciModeIdx] = max(dTTmatNorm(:,minBinIdx:end),[],2);
        saturatedSet = find(iciModeIdx <= 3);
        
        for iSat = 1:length(saturatedSet)
            if ~isnan(mean(dTTmatNorm(saturatedSet(iSat),1:end)));
                thisDTTsmoothDiff = diff(dTTmatNorm(saturatedSet(iSat),1:end));
                minIdx = find(thisDTTsmoothDiff(2:end)>0.02, 1,'first')+1;
                if isempty(minIdx)
                    minIdx = 2;
                end
                [mVal,mTemp] = max(dTTmatNorm(saturatedSet(iSat),minIdx:end),[],2);
                iciModeIdx(saturatedSet(iSat)) = minIdx+mTemp-1;
            end
        end
        iciMode = p.barInt(iciModeIdx+minBinIdx-1) + p.barInt(2)./2;
        
        % iterate over spectra in bin
        for iS = 1:size(specClickTf_norm_short,1)
            if nSpec{i0}(iS)>100
                %         iciDist = exp(-pdist2(iciMode(iS),iciCentroids,'euclidean'));
                %         specCorr = exp(-pdist2(specClickTf_norm_short(iS,:),...
                %                 specCentroids,'correlation'));
                %         comboDist = specCorr.*iciDist;
                
                % iterate over possible clusters
                comboDist = zeros(size(Tfinal,1),1);
                comboDistSortSet = {};
                comboDistStd = zeros(size(Tfinal,1),1);
                comboDistHist = [];
                thisSpec = (specClickTf_norm_short(iS,:)-min(specClickTf_norm_short(iS,:)))./...
                    max(specClickTf_norm_short(iS,:)-min(specClickTf_norm_short(iS,:)));
                
                
                specCorr = {};
                iciModeDist = {};
                for iTF = 1:size(Tfinal,1)
                    % compute spectral similarities
                    
                    specCorr{iTF} = exp(-pdist2(diff(thisSpec),...
                        diff(normTfinalSpec{iTF,1},1,2),'correlation'));
                    % compute modal ici similarities
                    %             iciModeDist{iTF}  = (exp(-pdist2(dTTmatNorm(iS,minICI:maxICI),...
                    %                 Tfinal{iTF,2}(:,minICI:maxICI),'correlation')));
                    iciModeDist{iTF} = exp(-pdist2(iciMode(iS),Tfinal{iTF,4}','euclidean')./p.barInt(maxICI));
                    % compute combined similarity scores
                    % [~,sortSpecScoreIdx] = sort(iciModeDist,'descend');
                    [comboDistSortSet{iTF},comboDistSortSetIdx{iTF}] = sort(specCorr{iTF}.*iciModeDist{iTF},'descend');
                    
                    %bestN = round(.5*length(comboDistSortSetIdx{iTF}));
                    %sumOfSq(iTF,1) = sum(comboDistSortSet{iTF}(1:bestN).^2)./bestN;%./length(comboDistSortSet{iTF});
                    % comboDistSortSet{iTF} = specCorr;
                    % comboDistHist(iTF,:) = histc(comboDistSortSet{iTF},0:.1:1);
                end
                pCut = prctile(cell2mat(comboDistSortSet),90);
                matchedSpec = [];
                nAboveCutoff = [];
                matchSet = [];
                prunedMatch = [];
                prunedMatchCount = [];
                for iPc = 1:size(comboDistSortSet,2)
                    keepSet = comboDistSortSet{iPc}(comboDistSortSet{iPc}>=pCut);
                    prunedMatch(iPc,1) = nansum(keepSet.^2)/length(keepSet);
                    prunedMatchCount(iPc,1) = length(keepSet);
                    % lowCount = find(prunedMatchCount(iPc,1)<(length(comboDistSortSet{iPc})*.1));
                    if prunedMatchCount(iPc,1)<(length(comboDistSortSet{iPc})*.1)
                        prunedMatch(iPc,1) = NaN;
                    end
                    matchSet{iPc} = 1:length(comboDistSortSet{iPc});
                    % nAboveCutoff(iPc) = length(matchSet{iPc});
                    comboDist(iPc) = nanmean(comboDistSortSet{iPc}(matchSet{iPc}));
                    matchedSpec(iPc,:) = nanmean(normTfinalSpec{iPc,1}(comboDistSortSetIdx{iPc}...
                        (1:min(length(comboDistSortSetIdx{iPc}),20)),:),1);
                    matchedICI(iPc,:) = nanmean(Tfinal{iPc,2}(comboDistSortSetIdx{iPc}...
                        (1:min(length(comboDistSortSetIdx{iPc}),20)),:),1);
                    meanspecScore(iPc,1) = nanmean(specCorr{iPc});
                    meaniciScore(iPc,1) = nanmean(iciModeDist{iPc});
                    
                end
                [C,I] = max(prunedMatch);
                [bestSpec,bestSpecIdx] = max(meanspecScore);
                [bestICI,bestICIIdx] = max(meaniciScore);
                %         if (bestSpecIdx == bestICIIdx) && bestICIIdx ~= I
                %             I = bestSpecIdx;
                %             disp(sprintf('iteration %d',i0))
                %             disp('Overriding match because best ICI and best spec agree on something different')
                %         end
                matchedSpec = specCentroids;
                matchedICI  = iciCentroids;
                
                
                
                if visualize || manualValidate
                    if  C >= cMax
                        cStr = 'b';
                    else
                        cStr = 'r';
                    end
                    figure(111);clf
                    subplot(1,2,1)
                    plot(f(p.stIdx:p.edIdx),thisSpec,cStr,'linewidth',2);
                    hold on
                    if C >= cMax
                        plot(f(p.stIdx:p.edIdx),matchedSpec(I,:),'k')
                        %title(sprintf('%d good matches',nAboveCutoff(I)));
                        
                    end
                    title(sprintf('best spec:%d, %0.2f \n best ici:%d, %0.2f',bestSpecIdx,...
                        bestSpec,bestICIIdx,bestICI))
                    ylim([0,1])
                    xlim([10,65])
                    hold off
                    
                    subplot(1,2,2)
                    bar(p.barInt,dTTmatNorm(iS,:),cStr)
                    hold on
                    if C >= cMax
                        plot(p.barInt(1:maxICI),matchedICI(I,:),'k','linewidth',2)
                    end
                    title(gca,sprintf('Type %d, Overall Score %0.2f; Bin %d of %d',I,C,i0,length(sumSpec)))
                    xlim([0,max(p.barInt)])
                    xlim([0,.4])
                    if C < cMax
                        1;
                    else
                        1;
                    end
                    1;
                end
                
                if ~isempty(intersect(myInt,i0)) %&& iS == 1
                    if manualValidate && C>=.05
                        if C > .25
                            tempAnswer = I;
                        else
                            tempAnswer = nTemplates+1;
                        end
                        manualAnswer = NaN;
                        while isnan(manualAnswer)
                            manualAnswer = input(sprintf('Please enter ID number (temp ans = %d):',tempAnswer));
                            if isempty(manualAnswer)
                                manualAnswer = tempAnswer;
                            elseif manualAnswer>nTemplates+1 || manualAnswer<1
                                manualAnswer = NaN; 
                            end
                        end
                        Iman = manualAnswer;
                    else
                        Iman = nTemplates+1; % if no good matches, assign a max IDx
                        % pause(2);
                        
                    end
                else
                    Iman = NaN;% pause(2);
                end
                autoID{i0} = [autoID{i0},I];
                manVAuto = [manVAuto;[I, Iman, C]];
                autoScore{i0} = [autoScore{i0},C];
                labeledCountsAll(i0,I) = labeledCountsAll(i0,I)+...
                    round((percSpec{i0}(iS)./sum(percSpec{i0}))*cInt(i0));
            else
                autoID{i0} = [autoID{i0},nTemplates+1];
                manVAuto = [manVAuto;[I, Iman, C]];
                autoScore{i0} = [autoScore{i0},NaN];
                labeledCountsAll(i0,nTemplates+1) = labeledCountsAll(i0,nTemplates+1)+...
                    round((percSpec{i0}(iS)./sum(percSpec{i0}))*cInt(i0));
            end
            if length(manVAuto)==65
                1;
            end
        end
        pSpecAll(i0,1) = percSpec(i0);
        cAll(i0,1) = cInt(i0,1);
    end
    if mod(i0,500)==0
        fprintf('Done with iteration %d of %d\n',i0, length(sumSpec))
    end
end

%%
figure(191);clf
nTypes = size(labeledCountsAll,2);
[ha, pos] = tight_subplot(nTypes, 1,.01,.05,[.1 .03]);

for iP = 1:nTypes
    axes(ha(iP))
    gt0 = find(labeledCountsAll(:,iP)>0);
    plot(tInt(gt0,1),labeledCountsAll(gt0,iP)./1000,'.')
    if iP<nTypes
        set(gca,'xtickLabel',[])
    else
        datetick('x','mmm-yy','keepLimits','keepTicks')
    end
    xlim([min(tInt(:,1)),max(tInt(:,1))])
    ylimHi = get(gca,'ylim');
    ylim([0,ylimHi(2)])
end
pY = mylabel(gcf,'Clicks per bin (1000s)','FontSize',13,'xoff',.1);
pX = mxlabel(gcf,'Date (mmm-yy)','FontSize',13,'xoff',0,'yoff',.05);
saveas(191,fullfile(savDir,sprintf('%s_autoMatch_timeseries.fig',siteStr)))

%%
figure(193);clf
tsWeekStart = floor(min(tInt(:,1)) - weekday(min(tInt(:,1)))+1);
[ha2, pos2] = tight_subplot(nTypes, 1,.01,.05,[.1 .03]);
tsWeekEnd = ceil(max(tInt(:,1)) - weekday(max(tInt(:,1)))+7);
tsWeekBins = tsWeekStart:7:tsWeekEnd;
[~,ItInt] = histc(tInt(:,1),tsWeekBins);
for iP = 1:nTypes
    axes(ha2(iP))
    cpw = zeros(length(tsWeekBins)-1,1);
    for iW = 1:length(tsWeekBins)-1
        cpw(iW,1) = sum(labeledCountsAll(ItInt==iW,iP)>0);
    end
    bar(tsWeekBins(1:end-1)+3.5,cpw)
    if iP<nTypes
        set(gca,'xtickLabel',[])
    else
        datetick('x','mmm-yy','keepLimits','keepTicks')
    end
    xlim([min(tInt(:,1)),max(tInt(:,1))])
    ylimHi = get(gca,'ylim');
    ylim([0,ylimHi(2)])
end
%%

manVAutoPrune = manVAuto(~isnan(manVAuto(:,2)),:);
manVAutoPrune = manVAutoPrune(manVAutoPrune(:,1)<nTemplates+1,:);
[C,I] = histc(manVAutoPrune(:,3),0:.1:1);
rightClass_perc = nan(10,1);
lenSet = [];
for iR = 1:length(C)
    autoSet = manVAutoPrune((I==iR),1);
    manSet = manVAutoPrune((I==iR),2);
    rightClass_perc(iR) = sum(all(autoSet==manSet,2))./length(autoSet);
    lenSet(iR) = length(autoSet);
end
save(fullfile(savDir,sprintf('%s_autoMatch',siteStr)),'siteStr',...
    'labeledCountsAll','pSpecAll','autoScore','nTemplates','Tfinal','typeFile',...
    'tInt','tsWeekBins','autoID','cAll','manualValidate','manVAuto','lenSet','rightClass_perc')
