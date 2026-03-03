function HierarchicalPriors_PLOT(parms)

% unpack parameter structure
v2struct(parms);

% start up Fieldtrip
addpath(FTdir);
ft_defaults

%% some parameters for plotting:
models2plot = 1:8;
latencytime = -maxlatency:1/fsNew:maxlatency;
tickTimes = [-1 -.5 0 .5 1];
tickIDs = dsearchn(latencytime',tickTimes');

cmap = brewermap([],'*RdBu');
lineColors = brewermap(7,'YlGnBu');% or 'RdPu'
lineColors(1,:) = [];
lineColors(1,:) = lineColors(1,:) - .1;

lineColorsLight = lineColors;
lineColorsDark = lineColors-0.2;
for c=1:numel(lineColors)
    if lineColorsLight(c) > 1
        lineColorsLight(c) = 1;
    end
    if lineColorsDark(c) < 0
        lineColorsDark(c) = 0;
    end
end

fs = 6;% fontsize

%% input and output folders
simstring = {'corr' ['pcaANDpcr_' num2str(nPCRcomps) 'comps']};
simstring = simstring{similarity+1};
indir = sprintf('%s%sresults%ssource_%s%sRSA%sstatistics%s%s_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,filesep,simstring,iternum,ceil(stimlen));

conditionnames = {'normal','invert','scram'};

load(fullfile(rootdir,'code','dRSA','ROIdefinitions'),'ROIdefinition');

ROIdefinition = ROIdefinition.names;
ROInames = vertcat(ROIdefinition{:});
ROInames = ROInames(ROI);

modellabels = {'pixelwise','OFmag','pixelwise motion','view-dependent body posture','view-invariant body posture','view-dependent body motion','view-invariant body motion','eye position'};

close(figure(1));
figure(1);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position',  [1 5 20 13]);

% initialize
peakMagnitudeROIs = nan(3,length(subjects),length(ROInames),length(modellabels)-1,2);
peakMagnitudeJKROIs = nan(3,length(subjects),length(ROInames),length(modellabels)-1,2);
peakLatencyJKROIs = nan(3,length(subjects),length(ROInames),length(modellabels)-1,2);
peakLatencyROIs = nan(3,length(ROInames),length(modellabels)-1,2);
peakAverageROIs = nan(3,length(subjects),length(ROInames),length(modellabels)-1,2);
dRSAGA = nan(length(subjects),3,length(ROInames),length(modellabels),length(latencytime));
peakLatencyTable = cell(18,8);
for icon = 1:3

    % initialize
    dRSASEM = zeros(length(ROInames),length(modellabels),length(latencytime));
    dRSAROIs = zeros(length(ROInames),length(modellabels),length(latencytime));
    signLineROIs = zeros(length(ROInames),length(modellabels),length(latencytime));
    for iROI = 1:length(ROInames)
        
        %% load data
        ROIname = ROInames{iROI};
        fn = sprintf('%s%cSTATS_%s_p%s_p%s_fisherz%d_%dHz_%s_dissimilarity%d_smMEG%d_smRDMneu%d_smRDMmod%d'...
            ,indir, filesep,  conditionnames{icon}, '001', '05', fisherz, fsNew, ROIname, dissimilarity, smoothNeuralData, smoothNeuralRDM, smoothModelRDM);
        load(fn,'dRSAall','signLine','peakMagnitude','peakAverage','peakLatencyJK','peakLatency','peakMagnitudeJK');%,'peakLOO');
        
        dRSASEM(iROI,:,:) = squeeze(stdErrMean(dRSAall,1));
        dRSAROIs(iROI,:,:) = squeeze(mean(dRSAall,'omitnan'));        
        signLineROIs(iROI,:,:) = signLine;
        
        if icon == 3
            peakMagnitudeROIs(icon,4:end,iROI,:,:) = peakMagnitude;
            peakMagnitudeJKROIs(icon,4:end,iROI,:,:) = peakMagnitudeJK;
            peakAverageROIs(icon,4:end,iROI,:,:) = peakAverage;
            peakLatencyJKROIs(icon,4:end,iROI,:,:) = peakLatencyJK;
            dRSAGA(4:end,icon,iROI,:,:) = dRSAall;
        else
            peakMagnitudeROIs(icon,:,iROI,:,:) = peakMagnitude;
            peakMagnitudeJKROIs(icon,:,iROI,:,:) = peakMagnitudeJK;
            peakAverageROIs(icon,:,iROI,:,:) = peakAverage;
            peakLatencyJKROIs(icon,:,iROI,:,:) = peakLatencyJK;
            dRSAGA(:,icon,iROI,:,:) = dRSAall;
        end
        peakLatencyROIs(icon,iROI,:,:) = peakLatency;        
        
    end% iROI loop
    
    % change zeros to NaN so it's easier to plot significance lines while leaving gaps for non-significant values:
    signLineROIs(signLineROIs==0) = NaN;
    ylims4lines = [.175 .065 .0176 .038 .038 .0176 .0176 .1];

    %% dRSA curves for motion models
    figure(1);
    
    modcount = 0;
    for imodel = 1:8
        modcount = modcount + 1;
        % averaged over video time
        subplot(3,8,modcount+(icon-1)*8);
        
        ylimit = ylims4lines(imodel);
        dRSA2plot = squeeze(dRSAROIs(:,imodel,:));
        SEM2plot = permute(repmat(squeeze(dRSASEM(:,imodel,:)),[1 1 2]),[2 3 1]);
        
        hold on;
        clear h
        
        boundedline(latencytime,dRSA2plot, SEM2plot , 'alpha','cmap', lineColors(1:6,:));
        
        for iROI = 1:size(dRSA2plot,1)
            
            h(iROI) = plot(latencytime,dRSA2plot(iROI,:),'color',lineColors(iROI,:),'lineWidth',1.5);
            
            plot(latencytime,(-ylimit*0.25-(ylimit/18*(iROI-1)))*squeeze(signLineROIs(iROI,imodel,:)),'lineWidth',3,'Color',lineColors(iROI,:));
            
        end
        
        plot([latencytime(1) latencytime(end)],[0 0],'--k');
        plot([0 0],[-1 1],'--k');
        set(gca,'xlim',[-1 1]);
        set(gca,'ylim',[-ylimit*.5 ylimit]);
        set(gca,'xtick',latencytime(tickIDs));
        set(gca,'xticklabel','');
        set(gca,'Fontsize',fs,'FontName','Helvetica');
        hold off

    end% iRDM

    % prepare table with peak latencies of subject average curves
    mod2plot = 0;
    for imodel = 1:length(models2plot)

        if imodel == 4
            dirs = 2;
        else
            mod2plot = mod2plot+1;
            dirs = 1;
        end

        for iROI = 1:length(ROInames)
            peakLatencyTable{iROI+(icon-1)*6,imodel} = num2str(1000*squeeze(peakLatencyROIs(icon,iROI,mod2plot,dirs)),'%.0f');
            if strcmp(peakLatencyTable{iROI+(icon-1)*6,imodel},'NaN')
               peakLatencyTable{iROI+(icon-1)*6,imodel} = '-';
            end
        end

    end
    
end% condition loop

% cd('XXX\HierarchicalPriors\figures');
% set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters HierarchicalPriors_ROIlagplots_all.eps

% writecell(peakLatencyTable,fullfile('XXX\HierarchicalPriors\results','peakLatencyTable.xls'));

%% STATS: normal vs inverted
% store peak magnitudes averaged in window 100 msec surrounding peak from Nat Commun study for comparision between normal and inverted
peakAverage4JASP = reshape(permute(squeeze(peakAverageROIs(1:2,:,1:4,[3 6 7],1)),[2 1 3 4]),length(subjects),[]);

% add average over ROIs
peakAverage4JASP(:,end+1:end+2) = squeeze(mean(peakAverageROIs(1:2,:,1:4,3,1),3))';
peakAverage4JASP(:,end+1:end+2) = squeeze(mean(peakAverageROIs(1:2,:,1:4,6,1),3))';
peakAverage4JASP(:,end+1:end+2) = squeeze(mean(peakAverageROIs(1:2,:,1:3,7,1),3))';

% store peak magnitudes
varnames = cell(2,4,3);
for icon = 1:2
    for iROI = 1:4
        modcount = 0;
        for imodel = [3 6 7]
            modcount = modcount + 1;
            varnames{icon,iROI,modcount} = [conditionnames{icon*2-1} '_' ROInames{iROI} '_' modellabels{imodel}];
        end
    end
end
varnames = reshape(varnames,1,[]);

% add ROI average
varnames{end+1} = 'norm_4ROI_OFdir';
varnames{end+1} = 'inv_4ROI_OFdir';
varnames{end+1} = 'norm_4ROI_depmot';
varnames{end+1} = 'inv_4ROI_depmot';
varnames{end+1} = 'norm_3ROI_invmot';
varnames{end+1} = 'inv_3ROI_invmot';

peakAverage4JASP = [varnames ; num2cell(peakAverage4JASP)];
writecell(varnames,fullfile('XXX\HierarchicalPriors\results','peakAverage4JASP.csv'));

%% STATS: normal vs scrambled preditive motion models
% % compare peak magnitudes at a priori peak latencies from Nat Commun
peakMagnitude4JASP = reshape(permute(squeeze(peakMagnitudeROIs([1 3],:,1:4,[3 6 7],1)),[2 1 3 4]),length(subjects),[]);

% add average over ROIs
peakMagnitude4JASP(:,end+1:end+2) = squeeze(mean(peakMagnitudeROIs([1 3],:,1:4,3,1),3))';
peakMagnitude4JASP(:,end+1:end+2) = squeeze(mean(peakMagnitudeROIs([1 3],:,1:4,6,1),3))';
peakMagnitude4JASP(:,end+1:end+2) = squeeze(mean(peakMagnitudeROIs([1 3],:,1:4,7,1),3))';

% remove first 3 participants because they miss scrambled condition
peakMagnitude4JASP(1:3,:) = [];
 
% create variable names
varnames = cell(2,4,3);
for icon = 1:2
    for iROI = 1:4
        modcount = 0;
        for imodel = [3 6 7]
            modcount = modcount+1;
            varnames{icon,iROI,modcount} = [conditionnames{icon*2-1} '_' ROInames{iROI} '_' modellabels{imodel}];
        end
    end
end
varnames = reshape(varnames,1,[]);

% add ROI average
varnames{end+1} = 'norm_4ROI_OFdir';
varnames{end+1} = 'scram_4ROI_OFdir';
varnames{end+1} = 'norm_4ROI_depmot';
varnames{end+1} = 'scram_4ROI_depmot';
varnames{end+1} = 'norm_4ROI_invmot';
varnames{end+1} = 'scram_4ROI_invmot';
 
peakMagnitude4JASP = [varnames ; num2cell(peakMagnitude4JASP)];
writecell(peakMagnitude4JASP,fullfile('XXX\HierarchicalPriors\results','peakMagnitude4JASP.csv'));

%% STATS: normal vs scrambled lagged poster models
% compare peak magnitudes at latencies of GA peaks
GA4peak = squeeze(mean(mean(mean(dRSAGA(4:40,[1 3],1:4,4:5,:),'omitnan'))));

% find peaks
[~,GApeakID] = max(GA4peak,[],2);

peakmagGA4JASP = zeros(37,2,4,2);
modcount = 0;
for imod = 4:5
    modcount = modcount+1;
    peakmagGA4JASP(:,:,:,modcount) = dRSAGA(4:40,[1 3],1:4,imod,GApeakID(modcount));
end

peakmagGA4JASP = reshape(peakmagGA4JASP,length(subjects)-3,[]);

% add average over ROIs
peakmagGA4JASP(:,end+1:end+2) = squeeze(mean(dRSAGA(4:40,[1 3],1:4,4,1),3));
peakmagGA4JASP(:,end+1:end+2) = squeeze(mean(dRSAGA(4:40,[1 3],1:4,5,1),3));
 
% create variable names
varnames = cell(2,4,2);
for icon = 1:2
    for iROI = 1:4
        modcount = 0;
        for imodel = 4:5
            modcount = modcount+1;
            varnames{icon,iROI,modcount} = [conditionnames{icon*2-1} '_' ROInames{iROI} '_' modellabels{imodel}];
        end
    end
end
varnames = reshape(varnames,1,[]);

% add ROI average
varnames{end+1} = 'norm_4ROI_deppos';
varnames{end+1} = 'scram_4ROI_deppos';
varnames{end+1} = 'norm_4ROI_invpos';
varnames{end+1} = 'scram_4ROI_invpos';

peakmagGA4JASP = [varnames ; num2cell(peakmagGA4JASP)];
writecell(peakmagGA4JASP,fullfile('XXX\HierarchicalPriors\results','peakmagGA_normVSscram_4JASP.csv'));

%% STATS: replication Nat Commun
% store peak latencies
peakLatencyJK4JASP_normal = reshape(squeeze(peakLatencyJKROIs(1,:,1:4,[3 6 7],1)),length(subjects),[]);

%% scatter plot figures of condition by condition
% normal vs inverted condition, predictive motion peaks
xlim = [min(peakAverage4JASP(:,25:2:end)) ; max(peakAverage4JASP(:,25:2:end))];
ylim = [min(peakAverage4JASP(:,26:2:end)) ; max(peakAverage4JASP(:,26:2:end))];
limxy = cat(3,xlim,ylim);
lim(1,:) = min(limxy(1,:,:),[],3);
lim(2,:) = max(limxy(2,:,:),[],3);
lim = lim*1.1;

close(figure(2))
figure(2);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position',  [1 5 8 10]);

for iplot = 1:3
    subplot(3,1,iplot)
    scatter(peakAverage4JASP(:,23+iplot*2),peakAverage4JASP(:,24+iplot*2),14,'LineWidth',.5,'MarkerEdgeColor','k','MarkerFaceColor',[.6 .6 .6]);
    hold on
    scatter(mean(peakAverage4JASP(:,23+iplot*2)),mean(peakAverage4JASP(:,24+iplot*2)),14,'LineWidth',.5,'MarkerEdgeColor','k','MarkerFaceColor','r');
    plot([-0.02 ; 0.05],[-0.02 ; 0.05],'k--','LineWidth',1);
    set(gca,'xlim',lim(:,iplot));
    set(gca,'ylim',lim(:,iplot));
    axis square
    box on
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',7,'FontName','Helvetica')
    set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0])
end

% cd('XXX\HierarchicalPriors\figures');
% set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters HierarchicalPriors_peakAvg_normVSinv_motion.eps

% normal vs scrambled condition, predictive motion peaks
xlim = [min(peakMagnitude4JASP(:,25:2:end)) ; max(peakMagnitude4JASP(:,25:2:end))];
ylim = [min(peakMagnitude4JASP(:,26:2:end)) ; max(peakMagnitude4JASP(:,26:2:end))];
limxy = cat(3,xlim,ylim);
lim(1,:) = min(limxy(1,:,:),[],3);
lim(2,:) = max(limxy(2,:,:),[],3);
lim = lim*1.1;

close(figure(3))
figure(3);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position',  [1 5 8 10]);

for iplot = 1:3
    subplot(3,1,iplot)
    scatter(peakMagnitude4JASP(:,23+iplot*2),peakMagnitude4JASP(:,24+iplot*2),14,'LineWidth',.5,'MarkerEdgeColor','k','MarkerFaceColor',[.6 .6 .6]);
    hold on
    scatter(mean(peakMagnitude4JASP(:,23+iplot*2)),mean(peakMagnitude4JASP(:,24+iplot*2)),14,'LineWidth',.5,'MarkerEdgeColor','k','MarkerFaceColor','r');
    plot([-0.02 ; 0.05],[-0.02 ; 0.05],'k--','LineWidth',1);
    set(gca,'xlim',lim(:,iplot));
    set(gca,'ylim',lim(:,iplot));
    axis square
    box on
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',7,'FontName','Helvetica')
    set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0])
end

% cd('XXX\HierarchicalPriors\figures');
% set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters HierarchicalPriors_peakMag_normVSscram_motion.eps

% normal vs scrambled condition, lagged posture peaks
xlim = [min(peakmagGA4JASP(:,17:2:end)) ; max(peakmagGA4JASP(:,17:2:end))];
ylim = [min(peakmagGA4JASP(:,18:2:end)) ; max(peakmagGA4JASP(:,18:2:end))];
limxy = cat(3,xlim,ylim);
clear lim
lim(1,:) = min(limxy(1,:,:),[],3);
lim(2,:) = max(limxy(2,:,:),[],3);
lim = lim*1.1;

close(figure(3))
figure(3);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position',  [1 5 8 10]);

for iplot = 1:2
    subplot(3,1,iplot)
    scatter(peakmagGA4JASP(:,15+iplot*2),peakmagGA4JASP(:,16+iplot*2),14,'LineWidth',.5,'MarkerEdgeColor','k','MarkerFaceColor',[.6 .6 .6]);
    hold on
    scatter(mean(peakmagGA4JASP(:,15+iplot*2)),mean(peakmagGA4JASP(:,16+iplot*2)),14,'LineWidth',.5,'MarkerEdgeColor','k','MarkerFaceColor','r');
    plot([-0.05 ; 0.05],[-0.05 ; 0.05],'k--','LineWidth',1);
    set(gca,'xlim',lim(:,iplot));
    set(gca,'ylim',lim(:,iplot));
    axis square
    box on
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',7,'FontName','Helvetica')
    set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0])
end

% cd('XXX\HierarchicalPriors\figures');
% set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters HierarchicalPriors_peakMagGA_normVSscram_posture.eps


end