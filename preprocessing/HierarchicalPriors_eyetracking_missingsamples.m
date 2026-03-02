function HierarchicalPriors_eyetracking_missingsamples(parms)

% small script to count amount of missing data in eyetracker signal to report in manuscript

% unpack parameter structure
v2struct(parms);

% set paths
indirEye = '\\XXX\HierarchicalPriors\data\MEG\ELpreprocessed';

%% loop over subjects
missedsamp_persub = zeros(length(subjects),1);
missedtrial_persub = zeros(length(subjects),1);
subcount = 0;
for isub = subjects

    subcount = subcount + 1;

    fn= sprintf('%s%cSUB%02d',indirEye, filesep, isub);
    load(fn,'missedSampleNum','missedTrialNum');

    % size of tested data structure to compute percentage of missing samples for manuscript: [84 trials * 6 signals * 5200 samples] = 2620800
    numrun = sum(~cellfun(@isempty, missedSampleNum(isub,:)));
    missedsamp_perrun = zeros(numrun,1);
    missedtrial_perrun = zeros(numrun,1);
    for irun = 1:numrun

        missedsamp_perrun(irun) = 100*sum(missedSampleNum{isub,irun})/2620800;
        missedtrial_perrun(irun) = 100*missedTrialNum(isub,irun)/84;

    end% run loop

    missedsamp_persub(subcount) = mean(missedsamp_perrun);
    missedtrial_persub(subcount) = mean(missedtrial_perrun);

end% subject loop

end
