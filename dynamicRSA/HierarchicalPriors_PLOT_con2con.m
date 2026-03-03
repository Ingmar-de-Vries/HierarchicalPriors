function HierarchicalPriors_PLOT_con2con(parms)

% unpack parameter structure
v2struct(parms);

% start up Fieldtrip
addpath(FTdir);
ft_defaults

%% some parameters for plotting:
latencytime = -maxlatency:1/fsNew:maxlatency;

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

fs = 10;% fontsize

%% input and output folders
simstring = {'corr_con2con' ,'corr_con2con_partout2'};
load(fullfile(rootdir,'code','dRSA','ROIdefinitions'),'ROIdefinition');

ROIdefinition = ROIdefinition.names;
ROInames = vertcat(ROIdefinition{:});
ROInames = ROInames(ROI);
ROInames(5) = {'ROI-avg'};

comb2load = 'normal2invert';

% initialize
dRSAROIs = zeros(2,length(ROI)+1,length(latencytime));
dRSASEM = dRSAROIs;
signLineROIs = dRSAROIs;
slopesL = zeros(2,length(subjects),length(ROI)+1);
slopesR = slopesL;
CI95ROIs1 = zeros(2,length(ROI)+1);
for isimi = 1:2% similarity loop

    indir = sprintf('%s%sresults%ssource_%s%sRSA%sstatistics%s%s_%diterations_%dsec',...
        rootdir,filesep,filesep,atlas2use,filesep,filesep,filesep,simstring{isimi},iternum,ceil(stimlen));

    %% load data
    fn = sprintf('%s%cSTATS_%s_p%s_p%s_fisherz%d_%dHz_dissimilarity%d_smMEG%d_smRDMneu%d'...
        ,indir, filesep,  comb2load, '001', '05', fisherz, fsNew, dissimilarity, smoothNeuralData, smoothNeuralRDM);
    load(fn,'dRSAperROI','signLine','CI95peak_1tail','CI95peak_2tail','slopeL','slopeR');

    dRSASEM(isimi,:,:) = squeeze(stdErrMean(dRSAperROI,1));
    dRSAROIs(isimi,:,:) = squeeze(mean(dRSAperROI,'omitnan'));
    signLineROIs(isimi,:,:) = signLine;
    slopesL(isimi,:,:) = abs(slopeL);
    slopesR(isimi,:,:) = abs(slopeR);

    CI95ROIs1(isimi,:) = CI95peak_1tail;
    CI95ROIs2(isimi,:,:) = CI95peak_2tail;
    [~,peakavg] = max(dRSAROIs,[],3);
    peakavg = latencytime(peakavg);

end% similarity loop

% temporally smooth curves for plotting
for isimi = 1:2
    dRSAROIs(isimi,:,:) = ft_preproc_smooth(squeeze(dRSAROIs(isimi,:,:)),9);
end% similarity loop

% change zeros to NaN so it's easier to plot significance lines while leaving gaps for non-significant values:
signLineROIs(signLineROIs==0) = NaN;

%% dRSA curves
close(figure(1));
figure(1);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position',  [1 5 16 10]);

for isimi = 1:2
    subplot(2,3,1+(isimi-1)*3)

    SEM2plot = repmat(squeeze(dRSASEM(isimi,5,:)),[1 2]);

    hold on;

    % plot average plus SEM shading
    boundedline(latencytime,squeeze(dRSAROIs(isimi,5,:)), SEM2plot , 'alpha','cmap', [0 0 0]);

    % plot average line again on top of boundedline but thicker
    plot(latencytime,squeeze(dRSAROIs(isimi,5,:)),'color',[0 0 0],'lineWidth',1.5);

    ylimit = max(max(squeeze(dRSAROIs(:,5,:)))) + max(max(squeeze(dRSASEM(:,5,:))));

    % plot statistics bars
    plot(latencytime,-ylimit*0.05*squeeze(signLineROIs(isimi,5,:)),'lineWidth',3,'Color',[0 0 0]);

    plot([latencytime(1) latencytime(end)],[0 0],'--k');
    plot([0 0],[-1 1],'--k');
    set(gca,'xlim',[-1 1]);
    set(gca,'ylim',[-ylimit*0.1 ylimit]);
    tickTimes = [-1 -.5 0 .5 1];
    tickIDs = dsearchn(latencytime',tickTimes');
    set(gca,'xtick',latencytime(tickIDs));
    set(gca,'Fontsize',8,'FontName','Helvetica');
    hold off

    if isimi == 1
        ylabel('dRSA [corr]','Fontsize',fs,'FontName','Helvetica');
    elseif isimi == 2
        ylabel('dRSA [part corr]','Fontsize',fs,'FontName','Helvetica');
    end

    % plot per ROI lines
    subplot(2,3,2+(isimi-1)*3)
    dRSAmax = max(squeeze(dRSAROIs(isimi,:,:)),[],2);

    hold on;

    for iROI = ROI
        plot(latencytime,squeeze(dRSAROIs(isimi,iROI,:))/dRSAmax(iROI),'color',lineColors(iROI,:),'lineWidth',1.5);
    end

    plot([latencytime(1) latencytime(end)],[0 0],'--k');
    plot([0 0],[-1 1],'--k');
    set(gca,'xlim',[-1 1]);%[-.05 .05]);
    set(gca,'ylim',[-0.1 1]);%[ylimit*0.3 ylimit]);
    tickTimes = [-1 -.5 0 .5 1];
    tickIDs = dsearchn(latencytime',tickTimes');
    set(gca,'xtick',latencytime(tickIDs));
    set(gca,'ytick',[0 .5 1]);
    set(gca,'YTickLabel',{'0','0.5','1.0'});
    set(gca,'Fontsize',8,'FontName','Helvetica');
    hold off

    % plot per ROI lines but zoomed in
    subplot(2,3,3+(isimi-1)*3)

    hold on;

    for iROI = ROI
        plot(latencytime,squeeze(dRSAROIs(isimi,iROI,:))/dRSAmax(iROI),'color',lineColors(iROI,:),'lineWidth',1.5);
    end

    plot([latencytime(1) latencytime(end)],[0 0],'--k');
    plot([0 0],[-1 1],'--k');
    set(gca,'xlim',[-.05 .05]);
    set(gca,'ylim',[.75 1]);
    tickTimes = [-0.05 0 .05];
    tickIDs = dsearchn(latencytime',tickTimes');
    set(gca,'xtick',latencytime(tickIDs));
    set(gca,'ytick',[0.8 0.9 1.0]);
    set(gca,'YTickLabel',{'0.8','0.9','1.0'});
    set(gca,'Fontsize',8,'FontName','Helvetica');
    hold off

end% similarity loop

% cd('XXX\HierarchicalPriors\figures');
% set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters HierarchicalPriors_normalVSinverted.eps

end