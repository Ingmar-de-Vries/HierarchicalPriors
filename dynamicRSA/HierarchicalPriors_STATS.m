function HierarchicalPriors_STATS(parms,~,iroi,icon,ithresh,~)

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

% subjects 1,2,3 miss the scrambled condition, so for those only run normal and inverted conditions
if icon == 3
    sub4stat(sub4stat<4) = [];
    subnum = length(sub4stat);
end

% divide ROIs
load(fullfile(rootdir,'code','dRSA','ROIdefinitions'),'ROIdefinition');
ROIname = ROIdefinition.names{iroi}{:};

% if relative dRSA, divide dRSA values by maximum possible dRSA value based on clean simulations without random noise
if dRSArelative == 1
    load(fullfile(indir,'..',['dRSAmax_' conditionnames{icon} '_' num2str(50) 'Hz']),'dRSAmax');
    dRSAmax = dRSAmax(models2test);
end

% for latency plots
latencytime = -maxlatency:1/fsNew:maxlatency;

% set thresholds for significance and create output filename
if ithresh == 1
    pthresh = 0.05;% threshold for cluster size test
    pthreshclust = 0.01;% single sample tests to determine which samples are included in cluster
elseif ithresh == 2
    pthresh = 0.05;% threshold for cluster size test
    pthreshclust = 0.001;% single sample tests to determine which samples are included in cluster
end

pthreshstring = num2str(pthresh);
pthreshstring(1:2) = [];
pthreshcluststring = num2str(pthreshclust);
pthreshcluststring(1:2) = [];

fnSTATS = sprintf('%s%cSTATS_%s_p%s_p%s_fisherz%d_%dHz_%s_dissimilarity%d_smMEG%d_smRDMneu%d_smRDMmod%d',...
    outdir, filesep, conditionnames{icon}, pthreshcluststring, pthreshstring, fisherz, fsNew, ROIname,...
    dissimilarity, smoothNeuralData, smoothNeuralRDM, smoothModelRDM);

%% load data, combine and average subjects
dRSAallSub = zeros(length(sub4stat),length(models2test),length(latencytime));
subcount = 0;
for isub = sub4stat
    subcount = subcount+1;

    fn = sprintf('%s%cSUB%02d_%s_%dHz_%s_dissimilarity%d_smMEG%d_smRDMneu%d_smRDMmod%d',...
        indir, filesep, isub, conditionnames{icon}, fsNew, ROIdefinition.names{iroi}{:},...
        dissimilarity, smoothNeuralData, smoothNeuralRDM, smoothModelRDM);
    load(fn,'dRSAall');

    % apply fisher transform if values are correlation
    if fisherz
        dRSAallSub(subcount,:,:) = atanh(dRSAall');
    else
        dRSAallSub(subcount,:,:) = dRSAall';
    end
end

dRSAall = dRSAallSub;
clear dRSAallSub

% if relative dRSA, divide dRSA values by maximum possible dRSA value based on clean simulations without random noise
if dRSArelative
    dRSAall = dRSAall ./ repmat(dRSAmax,[size(dRSAall,1) 1 size(dRSAall,3)]);
end

%% determine peak value and timing, using jackknife approach
MODwith2peaks = 3;% test for a positive and negative peak in OFdir

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
previouspeaks = squeeze(previouspeaks(:,:,iroi));% select current ROI

% The above peak latency indices are based on an fs of 100 Hz, so adjust if necessary
previouspeaks = round(previouspeaks/(100/fsNew));

% how wide is window to search for peak in samples
peakrange = round(0.05*fsNew);

% for OF dir limit search to positive or negative peak, i.e., don't let peak range overlap with the other peak. Cutoff is the trough between peaks
% in previous study. Only for first ROIs because last two ROIs don't show two significant peaks effects in previous and current study.
troughsOFdir = [109 109 107 105 101 101];
troughsOFdir = round(troughsOFdir(iroi)/(100/fsNew));% The above peak latency indices are based on an fs of 100 Hz, so adjust if necessary

% initialize
peakMagnitudeJK = zeros(subnum,length(models2test)-1,2);
peakMagnitude = zeros(subnum,length(models2test)-1,2);
peakAverage = zeros(subnum,length(models2test)-1,2);
peakLatencyJK = zeros(subnum,length(models2test)-1,2);
peakLatency = zeros(length(models2test)-1,2);
for imodel = 1:size(dRSAall,2)-1% skip eyetracker RDM

    if ~any(MODwith2peaks == imodel)% look for only single peak
        interval2test = [previouspeaks(1,imodel) - peakrange previouspeaks(1,imodel) + peakrange ; 1 1];
        dirs = 1;
    else% look for two peaks: a predictive and a lagged peak
        interval2test = [previouspeaks(1,imodel) - peakrange previouspeaks(1,imodel) + peakrange ; previouspeaks(2,imodel) - peakrange previouspeaks(2,imodel) + peakrange];
        dirs = 1:2;
    end

    % correct OF dir for 2 peaks
    if imodel == 3 && iroi < 5
        if interval2test(1,2) > troughsOFdir-1
            interval2test(1,2) = troughsOFdir-1;
        end
        if interval2test(2,1) < troughsOFdir+1
            interval2test(2,1) = troughsOFdir+1;
        end
    end

    % cut off interval if on edge
    interval2test(interval2test > length(latencytime)) = length(latencytime);
    interval2test(interval2test < 1) = 1;

    for idir = dirs% 1 = predictive peak, 2 = lagged peak

        for isub = 1:subnum

            curveSub = squeeze(dRSAall(isub,imodel,:));

            % average over interval
            peakAverage(isub,imodel,idir) = mean(curveSub(interval2test(idir,1):interval2test(idir,2)));

        end% isub loop

    end% idir loop

    % Alternatively, use a priori peaks and determine magnitude at exactly that that latency per subject
    % or non-jackknifed peak magnitude at a priori latencies from Nat Commun
    for idir = dirs% 1 = predictive peak, 2 = lagged peak
        peakID = previouspeaks(idir,imodel);
        peakMagnitude(:,imodel,idir) = dRSAall(:,imodel,peakID);
    end

    % for getting peak latency values from average peak, or using jackknifing to retrieve individual subject peak latencies:
    if ~any(MODwith2peaks == imodel)% look for only single peak
        interval2test = [1 length(latencytime) ; 1 1];
        dirs = 1;
    else% look for two peaks: a predictive and a lagged peak
        interval2test = [1 troughsOFdir ; troughsOFdir length(latencytime)];
        dirs = 1:2;
    end

    curveAvg = squeeze(mean(dRSAall(:,imodel,:)));
    for idir = dirs% 1 = predictive peak, 2 = lagged peak

        [peakVAL,peakID] = findpeaks(curveAvg(interval2test(idir,1):interval2test(idir,2)));

        if isempty(peakVAL)% if no peak found... select maximum
            [peakVAL,peakID] = max(curveAvg(interval2test(idir,1):interval2test(idir,2)));
        end

        % find highest of those
        [~,maxID] = max(peakVAL);
        peakID = peakID(maxID);
        peakID = peakID+interval2test(idir,1)-1;

        % store individual subject peak magnitudes
        peakLatency(imodel,idir) = latencytime(peakID);

        % and for jackknifed peak latencies and magnitudes:
        for isub = 1:subnum
            idx = 1:size(dRSAall,1) ~= isub;% leave-one-out
            curveJK = squeeze(mean(dRSAall(idx,imodel,:)));

            [peakVAL,peakID] = findpeaks(curveJK(interval2test(idir,1):interval2test(idir,2)));

            if isempty(peakVAL)% if no peak found... select maximum
                [peakVAL,peakID] = max(curveJK(interval2test(idir,1):interval2test(idir,2)));
            end

            % find highest of those
            [peakVAL,maxID] = max(peakVAL);
            peakID = peakID(maxID);
            peakID = peakID+interval2test(idir,1)-1;

            % store jackknifed peak latencies and dRSA values at those latencies
            peakMagnitudeJK(isub,imodel,idir) = peakVAL;
            peakLatencyJK(isub,imodel,idir) = latencytime(peakID);

        end% subject loop

    end% dir loop

    % Retrieve individual-subject latencies from subsample jackknife averages (Smulders 2010 Psychophysiology)
    % more of a 'how did this subject's peak influence the grand average'-measure, but that's the best we can get, and it allows for running normal stats
    n = size(peakLatencyJK,1);
    peakLatencyJK(:,imodel,:) = repmat(n .* mean(squeeze(peakLatencyJK(:,imodel,:))),n,1) - (n-1) .* squeeze(peakLatencyJK(:,imodel,:));
    peakMagnitudeJK(:,imodel,:) = repmat(n .* mean(squeeze(peakMagnitudeJK(:,imodel,:))),n,1) - (n-1) .* squeeze(peakMagnitudeJK(:,imodel,:));

end% iRDM loop

%% Statistics
chance = 0;% i.e. compare against a correlation of zero

%time series
testtype = 'multivar_timeseries';%univar_timeseries, univar_timefreq, univar_timefreqchan, multivar_timeseries, multivar_timefreq
signposLine = false(size(dRSAall,2),size(dRSAall,3));
signnegLine = false(size(dRSAall,2),size(dRSAall,3));

for imodel = 1:size(dRSAall,2)

    addinfo.time = latencytime;
    addinfo.tail = 1;% only positive tail
    [~,signposLine(imodel,:),signnegLine(imodel,:)] = HierarchicalPriors_statsFT(squeeze(dRSAall(:,imodel,:)),addinfo,[],testtype,pthresh,pthreshclust,chance);

end
signLine = logical(signposLine+signnegLine);

save(fnSTATS,'dRSAall','signLine','ROIname','peakMagnitude','peakMagnitudeJK','peakLatencyJK','peakLatency','peakAverage');

end