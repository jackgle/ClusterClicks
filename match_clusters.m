% match_clusters.m

% Takes output from cluster_bins.m and assigns click type labels to it
% using the output fom composite_clusters.m 

clearvars
siteStr = 'HZ';
savDir = 'H:\HAT02-03A\ClusterBins\Temp';
fileToMatch = 'H:\HAT02-03A\ClusterBins\HAT02and3A_cluster_bins.mat';
typeFile = 'H:\HAT02-03A\ClusterBins\Temp\HAT02-03_auto_95_typesHR.mat';

manualVerify = 0; % set to 1 if you want to manually assign an ID to bins with no good matches.
% else make it 0

% load data and templates
load(fileToMatch);
load(typeFile,'Tfinal');
% Note: If you don't want to match all the templates in Tfinal, prune it
% here, for example
% Tfinal = Tfinal([1,3,5],:);% only match templates 1, 3 and 5


% Some settings about what part of the spectra to compare.
p.edIdx = 121;
nTemplates = size(Tfinal,1);
cMax = 0.3;
visualize = 1;
minICI = 1;
maxICI = 31;
minNumClicks = 100; % don't try to classify bins with fewer than this number of clicks


% Check for output directory and make it if it doesn't exist
if ~exist(savDir,'dir')
    mkdir(savDir)
end

% Initialize variables
autoScore = cell(size(tInt,1),1);
autoID = cell(size(tInt,1),nTemplates);
labeledCountsAll = zeros(size(tInt,1),nTemplates+1);
pSpecAll = cell(size(tInt,1),1);
cAll = zeros(size(tInt,1),1);

%%
specCentroids = [];
iciCentroids = [];%nan(size(Tfinal,1),1);
nTf = ones(size(Tfinal,1),1);
normTfinalSpec = [];

for iTF = 1:size(Tfinal,1)
    iciCentroids(iTF,:) = nanmean(Tfinal{iTF,2});
    specCentroid1 = nanmean(Tfinal{iTF,1}(:,p.stIdx:p.edIdx));
    specCentroids(iTF,:) = (specCentroid1-min(specCentroid1))./max(specCentroid1-min(specCentroid1));
    nTf(iTF) = size(Tfinal{iTF},1);
    normTfinalSpec1 = Tfinal{iTF,1}(:,p.stIdx:p.edIdx);
    normTfinalSpec1 = normTfinalSpec1-repmat(min(normTfinalSpec1,[],2),1,length(p.stIdx:p.edIdx));
    normTfinalSpec{iTF,1} = normTfinalSpec1./repmat(max(normTfinalSpec1,[],2),1,length(p.stIdx:p.edIdx));
end


% if you are going to manually check, make a plot of all the click types to
% compare to.
if manualVerify
    figure(222);clf;
    subplot(1,2,1)
    plot(f(p.stIdx:p.edIdx),specCentroids)
    legend(num2str([1:nTemplates]'))
    xlim([10,65])
    subplot(1,2,2)
    plot(p.barInt(1:maxICI),iciCentroids)
end

minBinIdx = 1;
for i0 = 1:length(sumSpec)
    if cInt(i0)>1    
        % normalize mean spectra for this bin 
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
        
        % this statement handles saturated ICI distributions, where large
        % numbers of overlapping click trains produce a peak at low ICI.
        for iSat = 1:length(saturatedSet)
            if ~isnan(mean(dTTmatNorm(saturatedSet(iSat),1:end)))
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
            if nSpec{i0}(iS)>minNumClicks  % doesn't classify bins with fewer than X clicks
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
                    % Compute spectral and ICI similarities for each
                    % possible template.
                    
                    % There have been many variations of this, as evidenced
                    % by the commented code. Left in to provide alternate
                    % options.
                    
                    specCorr{iTF} = exp(-pdist2(diff(thisSpec),...
                        diff(normTfinalSpec{iTF,1},1,2),'correlation'));
                    % compute modal ici similarities
                    % iciModeDist{iTF}  = (exp(-pdist2(dTTmatNorm(iS,minICI:maxICI),...
                    % iTfinal{iTF,2}(:,minICI:maxICI),'correlation')));
                    iciModeDist{iTF} = exp(-pdist2(iciMode(iS),Tfinal{iTF,4}','euclidean')./p.barInt(maxICI));
                    % compute combined similarity scores
                    % [~,sortSpecScoreIdx] = sort(iciModeDist,'descend');
                    [comboDistSortSet{iTF},comboDistSortSetIdx{iTF}] = sort(specCorr{iTF}.*iciModeDist{iTF},'descend');
                    
                    % bestN = round(.5*length(comboDistSortSetIdx{iTF}));
                    % sumOfSq(iTF,1) = sum(comboDistSortSet{iTF}(1:bestN).^2)./bestN;%./length(comboDistSortSet{iTF});
                    % comboDistSortSet{iTF} = specCorr;
                    % comboDistHist(iTF,:) = histc(comboDistSortSet{iTF},0:.1:1);
                end
                
                % Prune out poor matches and pick the best overall matching
                % template. This step is finicky. I would like to see it
                % replaced with a neural net that would be trained on the
                % composite clusters output and then handed novel examples
                % from the full dataset.
                
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
                
                
                % Show the match, and if it's a bad match and the user
                % wants to manually verify, let them do so.
                if visualize || manualVerify
                    if  C >= cMax
                        % If the match is better than the threshold
                        % similariity cuttof, color blue, otherwise red
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
                    xlim([0,p.barInt(maxICI)])
                    drawnow
                    
%                     if C < cMax  % code for setting breakpoints during
%                     development
%                         1;
%                     else
%                         1;
%                     end
%                     1;
                end
                
                if C < cMax % if it's a bad match, have the user weigh in.
                    if manualVerify && C>=.05
                        if C > .25  
                            % If it's somewhat certain about the match, 
                            % it will suggest the best matching classification
                            tempAnswer = I;
                        else
                            % Otherwise it will sugges the unknown category.
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
                        I = manualAnswer;
                    else
                        I = nTemplates+1; % if no good matches, assign a max IDx
                        % pause(2);
                        
                    end
                else
                    % pause(2);
                end
                autoID{i0} = [autoID{i0},I];
                autoScore{i0} = [autoScore{i0},C];
                labeledCountsAll(i0,I) = labeledCountsAll(i0,I)+...
                    round((percSpec{i0}(iS)./sum(percSpec{i0}))*cInt(i0));
            else
                autoID{i0} = [autoID{i0},nTemplates+1];
                autoScore{i0} = [autoScore{i0},NaN];
                labeledCountsAll(i0,nTemplates+1) = labeledCountsAll(i0,nTemplates+1)+...
                    round((percSpec{i0}(iS)./sum(percSpec{i0}))*cInt(i0));
            end
        end
        pSpecAll(i0,1) = percSpec(i0);
        cAll(i0,1) = cInt(i0,1);
    end
    if mod(i0,500)==0
        fprintf('Done with iteration %d of %d\n',i0, length(sumSpec))
    end
end

%% Rough plotting to check the results
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
save(fullfile(savDir,sprintf('%s_autoMatch',siteStr)),'siteStr',...
    'labeledCountsAll','pSpecAll','autoScore','nTemplates','Tfinal','typeFile',...
    'tInt','tsWeekBins','autoID','cAll','manualVerify')


