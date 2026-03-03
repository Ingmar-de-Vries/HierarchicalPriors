function HierarchicalPriors_STATS_bootstrap(parms,~,~,~,~,~)

% unpack parameter structure
v2struct(parms);

% start up Fieldtrip
addpath(FTdir);
ft_defaults

%% input and output folders
simstring = {'corr' ['pcaANDpcr_' num2str(nPCRcomps) 'comps']};
simstring = simstring{similarity+1};
indir = sprintf('%s%sresults%ssource_%s%sRSA%s%s_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,simstring,iternum,ceil(stimlen));
outdir = sprintf('%s%sresults%ssource_%s%sRSA%sstatistics%s%s_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,filesep,simstring,iternum,ceil(stimlen));
if ~exist(outdir,'dir')
    mkdir(outdir);
end

conditionnames = {'normal','invert','scram'};

% divide ROIs
load(fullfile(rootdir,'code','dRSA','ROIdefinitions'),'ROIdefinition');
% ROIname = ROIdefinition.names{iroi}{:};

% for latency plots
latencytime = -maxlatency:1/fsNew:maxlatency;

%% load data, combine and average subjects
dRSAallSub = nan(length(sub4stat),length(models2test),length(conditionnames),4,length(latencytime));
dRSAmaxAll = zeros(length(conditionnames),length(models2test));
for iroi = 1:4
    for icon = 1:length(conditionnames)

        sub4stat = parms.sub4stat;
        subcount = 0;
        
        % subjects 1,2,3 miss the scrambled condition, so for those only run normal and inverted conditions
        if icon == 3
            sub4stat([1 2 3]) = [];
            subcount = 3;
        end

        for isub = sub4stat
            subcount = subcount+1;

            fn = sprintf('%s%cSUB%02d_%s_%dHz_%s_dissimilarity%d_smMEG%d_smRDMneu%d_smRDMmod%d',...
                indir, filesep, isub, conditionnames{icon}, fsNew, ROIdefinition.names{iroi}{:},...
                dissimilarity, smoothNeuralData, smoothNeuralRDM, smoothModelRDM);
            load(fn,'dRSAall');

            % apply fisher transform if values are correlation
            if fisherz
                dRSAallSub(subcount,:,icon,iroi,:) = atanh(dRSAall');
            else
                dRSAallSub(subcount,:,icon,iroi,:) = dRSAall';
            end
        end% subject loop

        % if relative dRSA, divide dRSA values by maximum possible dRSA value based on clean simulations without random noise
        if dRSArelative
            load(fullfile(indir,'..',['dRSAmax_' conditionnames{icon} '_' num2str(50) 'Hz']),'dRSAmax');
            dRSAmaxAll(icon,:) = dRSAmax(models2test);
        end

    end% condition loop
end% ROI loop
dRSAall = dRSAallSub;
clear dRSAallSub

% if relative dRSA, divide dRSA values by maximum possible dRSA value based on clean simulations without random noise
if dRSArelative
    for icon = 1:length(conditionnames)
        for imod = 1:length(models2test)
            dRSAall(:,imod,icon,:,:) = dRSAall(:,imod,icon,:,:) ./ dRSAmaxAll(icon,imod);
        end% model loop
    end% condition loop
end

%% determine peak latency using bootstrap approach
% MODwith2peaks = 3;% test for a positive and negative peak in OFdir

% determine windows for peak estimation a priori based on previous found peaks reported in Nat Comm
clear previouspeaks
previouspeaks(1,:,:) = [111 112 112 112 121 151 ;...
    108 108 108 108 109 108 ;...
    89 89 90 92 101 101 ;...
    116 118 121 121 201 1 ;...
    113 112 104 107 111 101 ;...
    85 83 82 81 69 83 ;...
    57 57 37 49 112 36];
previouspeaks(2,3,:) = [128 128 124 128 106 110];% lagged peaks for OF direction

% The above peak latency indices are based on an fs of 100 Hz, so adjust if necessary
previouspeaks = round(previouspeaks/(100/fsNew));

% how wide is window to search for peak in samples
peakrange = round(0.5*fsNew);

% for OF dir limit search to positive or negative peak, i.e., don't let peak range overlap with the other peak. Cutoff is the trough between peaks
% in previous study. Only for first ROIs because last two ROIs don't show two significant peaks effects in previous and current study.
troughsOFdir = [109 109 107 105 101 101];
troughsOFdir = round(troughsOFdir/(100/fsNew));% The above peak latency indices are based on an fs of 100 Hz, so adjust if necessary
troughsOFdir = round(mean(troughsOFdir(1:4)));% average trough over first 4 ROIs since we also average dRSA curves over those ROIs for OF direction

%% comparing Normal vs Inverted for OF direction, view-dependent body motion, and view-invariant body motion
% first average over ROIs that showed significant dRSA curves in at least 1 condition,
% which means the first 4 ROIs for OF dir and view-dep mot, and the first 3 ROIs for view-invar mot
dRSA4test = squeeze(mean(dRSAall(:,[3 6],1:2,:,:),4));% OF dir and view-dep mot
dRSA4test(:,3,:,:) = squeeze(mean(dRSAall(:,7,1:2,1:3,:),4));% view-invar mot

% determine a priori peak based on previous results
peak4test = squeeze(round(mean(previouspeaks(1,[3 6],1:4),3)));
peak4test(3) = squeeze(round(mean(previouspeaks(1,7,1:3),3)));

% initialize variables
peakboot = zeros(size(dRSA4test,2),nboot,size(dRSA4test,3));
comboot = peakboot;
peakdiffboot = zeros(size(dRSA4test,2),nboot);% initialize bootstrapped peak latency differences
comdiffboot = peakdiffboot;% initialize bootstrapped centre of mass differences
for imodel = 1:size(dRSA4test,2)% only motion models with predictive peaks
    
    interval2test = [peak4test(imodel) - peakrange peak4test(imodel) + peakrange];
    interval2test(interval2test<1) = 1;
    if imodel == 1% look for predictive beak only until through from previous study because that's where lagged peak starts 
        interval2test(2) = troughsOFdir;
    end

    for iboot = 1:nboot

        randsubID = randsample(size(dRSA4test,1),size(dRSA4test,1),true);

        peakspercon = zeros(length(randsubID),2);
        compercon = zeros(length(randsubID),2);
        subcount = 0;
        for isub = randsubID'
            subcount = subcount + 1;

            curve2test = squeeze(dRSA4test(isub,imodel,:,interval2test(1):interval2test(2)));

            for icon = 1:2

                [peakVAL,peakID] = findpeaks(curve2test(icon,:));

                if isempty(peakVAL)% if no peak found... select maximum
                    [peakVAL,peakID] = max(curve2test(icon,:));
                end

                % find highest of those
                [~,maxID] = max(peakVAL);
                peakID = peakID(maxID);

                peakspercon(subcount,icon) = peakID+interval2test(1)-1;

            end% condition loop

            % Or find the centre of mass (centroid) instead of the peak, which takes into account non-symmetry of the curve
            x = 1:size(curve2test,2);
           
            % ignore negative values because don't make sense, preferrably call them zero, but in case the entire curve is negative that doesn't work so in that case shift entire curve to positive 
            for icon = 1:2
                if all(curve2test(icon,:) < 0)
                    curve2test(icon,:) = curve2test(icon,:)-min(curve2test(icon,:),[],2);
                else
                    curve2test(icon,:) = max(curve2test(icon,:),0);
                end
            end% condition loop

            numerator = trapz(x.*curve2test,2);
            denominator = trapz(curve2test,2);

            compercon(subcount,:) = (numerator ./ denominator)+interval2test(1)-1;

        end% subject loop

        peakboot(imodel,iboot,:) = mean(latencytime(round(peakspercon)));
        peakdiffboot(imodel,iboot) = mean(-diff(peakspercon,[],2))/fsNew;
        comboot(imodel,iboot,:) = mean(latencytime(round(compercon)));
        comdiffboot(imodel,iboot) = mean(-diff(compercon,[],2))/fsNew;

    end% bootstrap loop

end% model loop

% select one-sided 95% confidence interval
peakdiffboot = sort(peakdiffboot,2);
comdiffboot = sort(comdiffboot,2);
CI95peak_2tail = [peakdiffboot(:,0.025*nboot) peakdiffboot(:,0.975*nboot)];% two-tailed
CI95peak_1tail = peakdiffboot(:,0.95*nboot);% one-tailed
peakdiffavg = mean(peakdiffboot,2);
CI95com_2tail = [comdiffboot(:,0.025*nboot) comdiffboot(:,0.975*nboot)];% two-tailed
CI95com_1tail = comdiffboot(:,0.95*nboot);% one-tailed
comdiffavg = mean(comdiffboot,2);

%% comparing Normal vs Scrambled for view-dependent and view-invariant body posture
% first average over ROIs that showed significant dRSA curves in at least 1 condition,
% which means the first 4 ROIs for OF dir and view-dep mot, and the first 3 ROIs for view-invar mot
dRSA4test = zeros(size(dRSAall,1)-3,2,2,size(dRSAall,5));
dRSA4test(:,1,:,:) = squeeze(mean(dRSAall(4:end,4,[1 3],:,:),4));% view-dep posture
dRSA4test(:,2,:,:) = squeeze(mean(dRSAall(4:end,5,[1 3],2:4,:),4));% view-invar posture

% determine a priori peak based on previous results
peak4test = squeeze(round(mean(previouspeaks(1,[3 6],1:4),3)));
peak4test(3) = squeeze(round(mean(previouspeaks(1,7,1:3),3)));

% initialize variables
peakboot = zeros(size(dRSA4test,2),nboot,size(dRSA4test,3));
comboot = peakboot;
peakdiffboot = zeros(size(dRSA4test,2),nboot);% initialize bootstrapped peak latency differences
comdiffboot = peakdiffboot;% initialize bootstrapped centre of mass differences
for imodel = 1:size(dRSA4test,2)% only posture models
    
    interval2test = [peak4test(imodel) - peakrange peak4test(imodel) + peakrange];
    interval2test(interval2test<1) = 1;
    interval2test(interval2test>size(dRSA4test,4)) = size(dRSA4test,4);

    for iboot = 1:nboot

        randsubID = randsample(size(dRSA4test,1),size(dRSA4test,1),true);

        peakspercon = zeros(length(randsubID),2);
        compercon = zeros(length(randsubID),2);
        subcount = 0;
        for isub = randsubID'
            subcount = subcount + 1;

            curve2test = squeeze(dRSA4test(isub,imodel,:,interval2test(1):interval2test(2)));

            for icon = 1:2

                [peakVAL,peakID] = findpeaks(curve2test(icon,:));

                if isempty(peakVAL)% if no peak found... select maximum
                    [peakVAL,peakID] = max(curve2test(icon,:));
                end

                % find highest of those
                [~,maxID] = max(peakVAL);
                peakID = peakID(maxID);

                peakspercon(subcount,icon) = peakID+interval2test(1)-1;

            end% condition loop

            % Or find the centre of mass (centroid) instead of the peak, which takes into account non-symmetry of the curve
            x = 1:size(curve2test,2);
            
            % ignore negative values because don't make sense, preferrably call them zero, but in case the entire curve is negative that doesn't work so in that case shift entire curve to positive 
            for icon = 1:2
                if all(curve2test(icon,:) < 0)
                    curve2test(icon,:) = curve2test(icon,:)-min(curve2test(icon,:),[],2);
                else
                    curve2test(icon,:) = max(curve2test(icon,:),0);
                end
            end% condition loop

            numerator = trapz(x.*curve2test,2);
            denominator = trapz(curve2test,2);

            compercon(subcount,:) = (numerator ./ denominator)+interval2test(1)-1;

        end% subject loop

        peakboot(imodel,iboot,:) = mean(latencytime(round(peakspercon)));
        peakdiffboot(imodel,iboot) = mean(-diff(peakspercon,[],2))/fsNew;
        comboot(imodel,iboot,:) = mean(latencytime(round(compercon)));
        comdiffboot(imodel,iboot) = mean(-diff(compercon,[],2))/fsNew;

    end% bootstrap loop

end% model loop

% select one-sided 95% confidence interval
peakdiffboot = sort(peakdiffboot,2);
comdiffboot = sort(comdiffboot,2);
CI95peak_2tail = [peakdiffboot(:,0.025*nboot) peakdiffboot(:,0.975*nboot)];% two-tailed
CI95peak_1tail = peakdiffboot(:,0.95*nboot);% one-tailed
peakdiffavg = mean(peakdiffboot,2);
CI95com_2tail = [comdiffboot(:,0.025*nboot) comdiffboot(:,0.975*nboot)];% two-tailed
CI95com_1tail = comdiffboot(:,0.95*nboot);% one-tailed
comdiffavg = mean(comdiffboot,2);

end