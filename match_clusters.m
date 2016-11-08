clearvars
siteStr = 'MC';
savDir = 'F:\GOM_clickTypePaper_detections\TPWS\ClusterOct2016';
inFile = 'F:\GOM_clickTypePaper_detections\TPWS\ClusterOct2016\MC03_disk14_clusters_PG0_PR95_MIN100_MOD0_noFP_typesHR.mat';
load(inFile);

nTemplates = size(Tfinal,1);
cMax = .65;
visualize = 0;

autoScore = cell(size(tInt,1),1);
autoID = cell(size(tInt,1),nTemplates);
labeledCountsAll = zeros(size(tInt,1),nTemplates+1);
pSpecAll = cell(size(tInt,1),1);
cAll = zeros(size(tInt,1),1);

%%
specCentroids = [];
iciCentroids = nan(size(Tfinal,1),1);
for iTF = 1:size(Tfinal,1)
    iciCentroids(iTF) = mean(Tfinal{iTF,4});
    specCentroids(iTF,:) = mean(Tfinal{iTF,3});
end

minBinIdx = 1;
for i1 = 1:length(sumSpec)
    % normalize mean spectra for this bin and compute diff
    minSSsection = min(sumSpec{i1},[],2);
    specClickTf_minNorm = (sumSpec{i1} - ...
        minSSsection(:,ones(1,size(sumSpec{i1},2))));
    maxSSsection = max(specClickTf_minNorm,[],2);
    specClickTf_norm = specClickTf_minNorm./maxSSsection(:,ones(1,size(specClickTf_minNorm,2)));
    specClickTf_norm_short = specClickTf_norm(:,stIdx:edIdx);
    specClickTf_norm_short = diff(specClickTf_norm_short,1,2);
    
    % determine modal ici(s)
    dTTmat = vertcat(dTT{i1});
    % normalize ICI distributions
    dTTmatNorm1 = dTTmat./repmat(sum(dTTmat,2),1,size(dTTmat,2));
    dTTmatNorm = dTTmatNorm1./repmat(max(dTTmatNorm1,[],2),1,size(dTTmat,2));
    
    [~,iciModeIdx] = max(dTTmatNorm(:,minBinIdx:end),[],2);
    iciMode = p.barInt(iciModeIdx+minBinIdx-1) + p.barInt(2)./2;
    
    % iterate over spectra in bin
    for iS = 1:size(specClickTf_norm_short,1)
        iciDist = exp(-pdist2(iciMode(iS),iciCentroids,'euclidean'));
        specCorr = exp(-pdist2(specClickTf_norm_short(iS,:),...
                specCentroids,'correlation'));
        comboDist = specCorr.*iciDist;

%         % iterate over possible clusters
%         comboDist = zeros(size(Tfinal,1),1);
%         
%         for iTF = 1:size(Tfinal,1)
%             % compute spectral similarities
%             specCorr = exp(-pdist2(specClickTf_norm_short(iS,:),...
%                 Tfinal{iTF,3},'correlation'));
%             % compute modal ici similarities
%             %iciModeDist = exp(-pdist2(dTTmatNorm(iS,1:maxICI),Tfinal{iTF,2},'euclidean'));
%             iciModeDist = exp(-pdist2(iciMode(iS),Tfinal{iTF,4}','euclidean'));
%             % compute combined similarity scores
%             comboDistSort = specCorr.*iciModeDist;
%             comboDist(iTF) = nanmean(comboDistSort);
%         end
        [C,I] = max(comboDist);
        
        if C < cMax
            I = nTemplates+1; % if no good matches, assign a max IDx
        end
        autoID{i1} = [autoID{i1},I];
        autoScore{i1} = [autoScore{i1},C];
        labeledCountsAll(i1,I) = labeledCountsAll(i1,I)+...
            round((percSpec{i1}(iS)./sum(percSpec{i1}))*cInt(i1));
        
        if visualize
            if  C >= cMax
                cStr = 'b';
            else
                cStr = 'r';
            end
            figure(111);clf
            subplot(1,2,1)
            plot(f,specClickTf_norm(iS,:),cStr,'linewidth',2);
            hold on
            if C >= cMax
                plot(f,Tfinal{I,5},'k')
            end
            ylim([0,1])
            xlim([10,65])
            hold off
            
            subplot(1,2,2)
            bar(p.barInt,dTTmatNorm(iS,:),cStr)
            hold on
            if C >= cMax
                plot(p.barInt(1:maxICI),mean(Tfinal{I,2}),'k','linewidth',2)
            end
            title(gca,sprintf('Overall Score %0.2f',C))
            xlim([0,max(p.barInt)])
            if C < cMax
                1;
            end
        end
    end
    pSpecAll(i1,1) = percSpec(i1);
    cAll(i1,1) = cInt(i1,1);
end



save(fullfile(savDir,sprintf('All%s_autoMatch',siteStr)),'siteStr',...
    'labeledCountsAll','pSpecAll','autoScore','nTemplates','Tfinal')


