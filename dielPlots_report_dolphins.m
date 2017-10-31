% Script to make MPL report plots from dolphin classifier output

% REQUIRES: Tethys plotting code including visPresence and visWeeklyPlot

% I usually save a version of this with each report I do, so that I
% remember what configuration I ran, what my final file choice was, etc.
clearvars
close all

% Database is used to get nighttime for deployment
q = dbInit('Server','bandolero.ucsd.edu');
dbSpeciesFmt('Input', 'Abbrev', 'NOAA.NMFS.v1') 	% enable our species abbreviations
saveTF = 1;
load('E:\Data\John Reports\JAX12D\JAX12D\dolphins_diff\matchOutput\JAX_test2_autoMatch.mat')
load('E:\Data\John Reports\JAX12D\JAX12D\dolphins_diff\JAX_auto_95_typesHR.mat','tInt')

savePath = 'E:\Data\John Reports\JAX12D\JAX12D';
load(fullfile('E:\Data\John Reports\JAX12D\JAX12D','plotbase_JAX12D.mat'))

% Note, this is currently hard-coded so that you can force it to match 
% plots generated for other species. Could grab this info from Tethys
% though...
effort = [datenum([2015,7,2,0,0,0]),datenum([2015,11,4,0,0,0])];
effortRND = [datenum([2015,7,2,0,0,0]),datenum([2015,11,4,0,0,0])];


labCountM = zeros(size(labeledCountsAll,1),3);
% for m1=1:size(labeledCountsAll,1)
%     lT = sum(labeledCountsAll(m1,[1,3,4,6])>0);% Risso's
%     if lT>0
%         labCountM(m1,1) = 1;
%     end
% end

% Clusters are never perfect. If you have multiple clicks types that are
% duplicates, merge them here. For example, there might be multiple flavors
% of Risso's clicks. Also, leave out things you don't want matched, or put
% them all in an unknown category.
labCountM(:,1) = sum(labeledCountsAll(:,[1,5,6]),2)>0;% Risso's
labCountM(:,2)= labeledCountsAll(:,7)>0;% CTJ1
labCountM(:,3)= labeledCountsAll(:,8)>0;% CTJ2
labCountM(:,4)= sum(labeledCountsAll(:,[2,3,4]),2)>0;% unk


lAll = sum(labCountM,2);
% labCountM(:,5)=lAll>0;
pLat = 30+(09.036/60);
pLon = 360-79+(46.203/60);
pNightTime = dbDiel(q,pLat,pLon,effortRND(1),effortRND(2));
% illu = dbGetLunarIllumination(q, pLat, pLon, effort(1,1), effort(2,2), 30);
xtick = 3; %default 3
ytick = 14; %default 7
resolution=10; %hourly resolution
    

%remove small breaks (<54 min) in effort...
% this section is taken from Tethys' plotting tool, to mimic the same
% effort treatment, but it doesn't actually have any effect as
% currently implemented, because of the manually input effort above.
if ~isempty(effort)
    effort = sortrows(effort,1);
    total = size(effort,1);
    for j = total-1:-1:1
        f1 = effort(j,1);
        f2 = effort(j+1,1);
        s1 = effort(j,2);
        s2 = effort(j+1,2);
        isOverlap = (f2 <= s1 + 0.0375);
        if isOverlap
            if s1 < s2
                effort(j,2) = effort(j+1,2);
            end
            effort(j+1,:) = [];
        end
    end
end
pEffort = effort;

%% Iterate over species and make plots
for n=1:size(labCountM,2)
    
    species = num2str(n);
    call = 'clicks';
    subtype = 'delph';
    deployments = 'JAX12D';
    title = sprintf('Plot %d',n);
    rect = [0, 0, 500, 600];
    figh = figure('NumberTitle','off','Name',title,'Position',rect,...
        'Units','Pixel');
    % add diel information
    nightH = visPresence(pNightTime, 'Color', 'black', ...
        'LineStyle', 'none', 'Transparency', .15, ...
        'Resolution_m', 1/60,'DateTickInterval',ytick,...
        'DateRange',effortRND, 'HourTickInterval',...
        xtick);
    %lunarH = visLunarIllumination(illu);
    
    % make diel plot
    pDetections = [tInt(labCountM(:,n)==1)']';
    sightH = visPresence(pDetections(:,1), ...
        'Resolution_m', 5,'DateTickInterval',ytick,...
        'LineStyle', 'none', 'DateRange', effortRND,...
        'Effort', pEffort, 'Label', title,...
        'HourTickInterval', xtick);
    set(gca, 'YDir', 'reverse');  %upside down plot
    
    if saveTF
        % save image to file(both png and fig versions)
        fname_png = sprintf('%s\\%s_%s_%s_%s_diel.png',savePath,deployments,species,call,subtype);
        fname_fig = sprintf('%s\\%s_%s_%s_%s_diel.fig',savePath,deployments,species,call,subtype);
        
        set(figh,'PaperPositionMode','auto')
        print(figh,'-dpng',fname_png,'-r300')
        saveas(figh,fname_fig)
    end
    
    
    % make timeseries plot
    for iD = 1:length(plotbase.week)-1
        plotbase.cum_hrs(iD,1) = sum(tInt(labCountM(:,n)==1)>=plotbase.week(iD)...
            & tInt(labCountM(:,n)==1)<plotbase.week(iD+1))./12;
    end
    plotExtents = [plotbase.min_date;plotbase.max_date];
    plotbase.maxY = max(plotbase.cum_hrs);
    oblongR = [0, 0, 800, 300];
    figh2 = figure('NumberTitle','off','Name',title,'Position',oblongR,...
        'Units','Pixel');
    visWeeklyPlot(plotbase,plotExtents,1)
    if saveTF
        % save image to file (both png and fig versions)
        fname2_png = sprintf('%s\\%s_%s_%s_%s_effort.png',savePath,deployments,species,call,subtype);
        fname2_fig = sprintf('%s\\%s_%s_%s_%s_effort.fig',savePath,deployments,species,call,subtype);
        set(figh2,'PaperPositionMode','auto')
        print(figh2,'-dpng',fname2_png,'-r300')
        saveas(figh2,fname2_fig)
        
    end
end


