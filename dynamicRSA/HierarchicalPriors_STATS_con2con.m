function HierarchicalPriors_STATS_con2con(parms,~,~,~,ithresh,~)

% unpack parameter structure
v2struct(parms);

% start up Fieldtrip
addpath(FTdir);
ft_defaults

%% input and output folders
simstring = {'corr_con2con' ['pcr_con2con_' num2str(nPCRcomps) 'comps'],'corr_con2con_partout2'};
simstring = simstring{similarity+1};
indir = sprintf('%s%sresults%ssource_%s%sRSA%s%s_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,simstring,iternum,ceil(stimlen));
outdir = sprintf('%s%sresults%ssource_%s%sRSA%sstatistics%s%s_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,filesep,simstring,iternum,ceil(stimlen));
if ~exist(outdir,'dir')
    mkdir(outdir);
end

conditionnames = {'normal','invert','scram'};

% % divide ROIs
load(fullfile(rootdir,'code','dRSA','ROIdefinitions'),'ROIdefinition');

% for latency plots
latencytime = -maxlatency:1/fsNew:maxlatency;

%% load data, combine and average subjects
dRSAperROI = zeros(length(sub4stat),length(ROI),length(latencytime));
subcount = 0;
for isub = sub4stat
    subcount = subcount+1;

    for iroi = ROI

        fn = sprintf('%s%cSUB%02d_%s2%s_%dHz_%s_dissimilarity%d_smMEG%d_smRDMneu%d',...
            indir, filesep, isub, conditionnames{con2con(1)}, conditionnames{con2con(2)}, fsNew, ROIdefinition.names{iroi}{:},...
            dissimilarity, smoothNeuralData, smoothNeuralRDM);
        load(fn,'dRSAall');

        if fisherz
            dRSAall = atanh(dRSAall);
        end

        dRSAperROI(subcount,iroi,:) = dRSAall';%smooth(dRSAall',5);% smooth for improved peak estimation

    end

end

% add ROI-average
dRSAperROI(:,5,:) = squeeze(mean(dRSAperROI,2));

%% bootstrapping peak latency against zero
latencies = zeros(nboot,size(dRSAperROI,2));
coms = latencies;
interval2test = [-.05 .05];
interval2test = dsearchn(latencytime',interval2test');
interval2test = interval2test(1):interval2test(2);
time4test = latencytime(interval2test);
for iroi = 1:size(dRSAperROI,2)
    for iboot = 1:nboot

        randsubID = randsample(length(sub4stat),length(sub4stat),true);

        latency = zeros(length(randsubID),1);
        com = latency;
        subcount = 0;
        for isub = randsubID'
            subcount = subcount + 1;

            curve2test = squeeze(dRSAperROI(isub,iroi,interval2test));

            [peakVAL,peakID] = findpeaks(curve2test);

            if isempty(peakVAL)% if no peak found... select maximum
                [peakVAL,peakID] = max(curve2test);
            end

            % find highest of those
            [~,maxID] = max(peakVAL);
            peakID = peakID(maxID);

            latency(subcount) = time4test(peakID);

            % Or find the centre of mass (centroid) instead of the peak, which takes into account non-symmetry of the curve
            x = 1:length(curve2test);

            % ignore negative values because don't make sense, preferrably call them zero, but in case the entire curve is negative that doesn't work so in that case shift entire curve to positive
            if all(curve2test < 0)
                curve2test = curve2test-min(curve2test);
            else
                curve2test = max(curve2test,0);
            end

            numerator = trapz(x.*curve2test',2);
            denominator = trapz(curve2test',2);

            com(subcount) = time4test(round(numerator / denominator));

        end

        latencies(iboot,iroi) = mean(latency);
        coms(iboot,iroi) = mean(com);

    end% bootstrap loop
end% ROI loop

% select one-sided 95% confidence interval
latencies = sort(latencies);
coms = sort(coms);
CI95peak_2tail = [latencies(0.025*nboot,:) ; latencies(0.975*nboot,:)];
CI95peak_1tail = latencies(0.95*nboot,:);
CI95com_2tail = [coms(0.025*nboot,:) ; coms(0.975*nboot,:)];
CI95com_1tail = coms(0.95*nboot,:);
peakBST = mean(latencies);
comBST = mean(coms);

%% Statistics
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

fnSTATS = sprintf('%s%cSTATS_%s2%s_p%s_p%s_fisherz%d_%dHz_dissimilarity%d_smMEG%d_smRDMneu%d',...
    outdir, filesep,  conditionnames{con2con(1)}, conditionnames{con2con(2)}, pthreshcluststring, pthreshstring, fisherz, fsNew,...
    dissimilarity, smoothNeuralData, smoothNeuralRDM);

chance = 0;% i.e. compare against a correlation of zero

%time series
testtype = 'multivar_timeseries';%univar_timeseries, univar_timefreq, univar_timefreqchan, multivar_timeseries, multivar_timefreq

addinfo.time = latencytime;
addinfo.tail = 1;% only positive tail

signLine = zeros(size(dRSAperROI,2),length(latencytime));
for iroi = 1:size(dRSAperROI,2)
    [~,signposLine,signnegLine] = HierarchicalPriors_statsFT(squeeze(dRSAperROI(:,iroi,:)),addinfo,[],testtype,pthresh,pthreshclust,chance);

    signLine(iroi,:) = logical(signposLine+signnegLine);
end

save(fnSTATS,'dRSAperROI','signLine','CI95peak_2tail','CI95peak_1tail','CI95com_2tail','CI95com_1tail','peakBST','comBST');

end