function HierarchicalPriors_PLOT_simulations(parms)

% unpack parameter structure
v2struct(parms);

% start up Fieldtrip
addpath(FTdir);
ft_defaults

%% some parameters for plotting:
latencytime = -maxlatency:1/fsNew:maxlatency;
tickTimes = [-1 -.5 0 .5 1];
tickIDs = dsearchn(latencytime',tickTimes');

colors = brewermap(9,'*RdYlBu');
colors(1:5,:) = colors(1:5,:) - .1;% first 5 colors are a bit too bright
colors(10,:) = [0 0 0];% last color is eye position which can be black
fs = 6;% fontsize

%% input and output folders
simstring = {'corr' ['pcaANDpcr_' num2str(nPCRcomps) 'comps']};
simstring = simstring{similarity+1};
indir = sprintf('%s%sresults%ssource_%s%sRSA%ssimulations%s%s_%dlag_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,filesep,simstring,lag*1000,iternum*length(subjects),ceil(stimlen));

conditionnames = {'normal','invert','scram'};

%% load data
% simulated data with separately all single model implanted
dRSA = zeros(length(subjects),length(models2sim),length(latencytime),length(modelnames));
regborder = zeros(length(subjects),length(modelnames));
for isub = 1:length(subjects)
    
    % I forgot to store the proper scramble indices for the first 3 subjects so we can't use the scramble condition, but the other 2 conditions are usable
    if subjects(isub) < 4 && con2plot == 3
        continue
    end
    
    for implant = 1:length(models2sim)
        
        fn2load = sprintf('%s%cSUB%02d_%s_implant%d_%dHz_dissimilarity%d_smMEG%d_smRDMneu%d_smRDMmod%d',...
            indir, filesep, subjects(isub), conditionnames{con2plot}, models2sim(implant), fsNew, dissimilarity, smoothNeuralData, smoothNeuralRDM, smoothModelRDM);

        load(fn2load,'dRSAall');
        dRSA(isub,implant,:,:) = dRSAall;
        
    end
    
    zeroID = dsearchn(latencytime',0);
    for imodel = 1:size(dRSA,4)
        
        peak = max(dRSA(isub,imodel,:,imodel));
        
        % based on 32% of peak of corr = 1 (i.e., 10% shared variance)
        Lslope = find(flip(squeeze(dRSA(isub,imodel,1:zeroID,imodel))) < .3162*peak,1);
        Rslope = find(squeeze(dRSA(isub,imodel,zeroID:end,imodel)) < .3162*peak,1);
        if isempty([Lslope Rslope])
            avg = maxlatency*fsNew;
        else
            avg = ceil(mean([Lslope Rslope],'omitnan'));
        end
        regborder(isub,imodel) = avg;
        
    end
        
end% subject loop

% average over subjects, skipping the first three for the scram condition
if con2plot == 3
    dRSA(1:3,:,:,:) = [];
end

dRSA = squeeze(mean(dRSA));

% store regression borders for attenuating model autocorrelation - only if similarity = 0 (correlation)
if similarity == 0

    temp.subinvarmods = round(mean(regborder(:,1:9)));
    temp.subvarmods = squeeze(regborder(:,10));
    regborder = temp;
    save(fullfile(indir,'..','..',['regressionBorder_smRDM30msec_' num2str(fsNew) 'Hz_' conditionnames{con2plot}]),'regborder');

elseif similarity == 1

    % store max peak values per model for normalizing real results to maximum possible dRSA value - only if similarity = 1 (regression)
    dRSAmax = zeros(1,length(modelnames));
    for imodel = 1:size(dRSA,3)
        dRSAmax(imodel) = max(squeeze(dRSA(imodel,:,imodel)));
    end
    save(fullfile(indir,'..','..',['dRSAmax_' conditionnames{con2plot} '_' num2str(fsNew) 'Hz']),'dRSAmax');

end

%% lag line plots (averaged over video time), ROIs combined in 1 figure
close(figure(1));

figure(1);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 5 20 3.5]);

for imodel = 1:length(models2test)
    subplot(1,8,imodel);
    
    if similarity == 0
        ylimit = 1;
    else
        ylimit = 1.01*max(max(max(squeeze(dRSA))));
    end
    h = [];
    
    hold on;
    h(1) = plot([latencytime(1) latencytime(end)],[0 0],'--k');
    for implant = 1:size(dRSA,1)-1
        h(implant+1) = plot(latencytime,squeeze(dRSA(implant,:,models2test(imodel))),'Color',colors(implant,:),'lineWidth',1);
    end
    h(end+1) = plot(latencytime,squeeze(dRSA(end,:,models2test(imodel))),'--','Color',[.2 .2 .2],'lineWidth',1);
        
    plot([0 0],[-1 1.5],'--k');
    plot([lag lag],[-1 1],'--k');
    set(gca,'xlim',[-1 1]);
    set(gca,'ylim',[-ylimit/10 ylimit]);
    set(gca,'xtick',latencytime(tickIDs));
    set(gca,'xtickLabel','');
    set(gca,'ytick',[0 .25 .5 .75 1]);
    set(gca,'ytickLabel','');
    set(gca,'Fontsize',fs,'FontName','Helvetica');
    hold off
    
end% iRDM

% cd('G:\My Drive\Active projects\HierarchicalPriors\figures');
% % set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters HierarchicalPriors_simulations_CORRnormal.eps
% print -depsc -tiff -r600 -painters HierarchicalPriors_simulations_PCRnormal.eps
% print -depsc -tiff -r600 -painters HierarchicalPriors_simulations_CORRscram.eps
% print -depsc -tiff -r600 -painters HierarchicalPriors_simulations_PCRscram.eps

end