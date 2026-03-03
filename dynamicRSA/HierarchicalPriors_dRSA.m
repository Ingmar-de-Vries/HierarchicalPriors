function HierarchicalPriors_dRSA(parms,isub,iROI,~,~,~)

% unpack parameter structure
v2struct(parms);

% turn off warmings on cluster because we don't see them anyway and might take time to print
if runcluster
    warning('off')
end

% start up Fieldtrip
addpath(FTdir);
ft_defaults

%% input and output folders
indirMEG = fullfile(rootdir,'data','MEG','PreprocessedSensor');
indirMOD = fullfile(rootdir,'data','modelRDMs');
indirEYE = fullfile(rootdir,'data','MEG','ELpreprocessed');
indirPTB = fullfile(rootdir,'data','MEG','PTBoutput');
indirSource = fullfile(rootdir,'data','MEG','brainstorm_database','HierarchicalPriors','data',subjectKEY{isub,2},subjectKEY{isub,2});
indirAtlas = fullfile(rootdir,'data','MEG','brainstorm_database','HierarchicalPriors','anat',subjectKEY{isub,2});

% create directory for results
simstring = {'corr' ['pcaANDpcr_' num2str(nPCRcomps) 'comps']};
simstring = simstring{similarity+1};
outdir = sprintf('%s%sresults%ssource_%s%sRSA%s%s_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,simstring,iternum,ceil(stimlen));
if ~exist(outdir,'dir')
    mkdir(outdir);
end

load(fullfile(rootdir,'code','dRSA','ROIdefinitions'),'ROIdefinition');
conditionnames = {'normal','invert','scram'};

%% load and prepare neural data
fn2load = sprintf('%s%cpreprocessedSensor_SUB%02d_%s.mat',indirMEG,filesep,isub,subjectKEY{isub,2});
load(fn2load,'FTdata'); % load in ft format

dataMEG = FTdata; clear FTdata
BSTevents = dataMEG.BSTevents;

% smooth with sliding window
if smoothNeuralData>0
    fsample = round(1/(dataMEG.time(2)-dataMEG.time(1)));
    smoothingSamples = smoothNeuralData / (1000/fsample);%transform smoothing in msec to smoothing in samples, because that's what ft_preproc_smooth uses
    for i = 1:size(dataMEG.trial,1)
        dataMEG.trial(i,:,:) = ft_preproc_smooth(squeeze(dataMEG.trial(i,:,:)), smoothingSamples);
    end
end

% downsample
if ~any(fsNew==[0 500])
    cfg = [];
    cfg.resamplefs = fsNew;
    cfg.detrend = 'no';
    dataMEG = ft_resampledata(cfg, dataMEG);
else
    fsNew = round(1/(dataMEG.time(2)-dataMEG.time(1)));
end

% For catch trials, find test onset and change all data after test onset to NaNs,
% this way, when averaging over trials using mean with the omitnan flag,
% only the pre-test interval is used for those trials
% load PTB output with test onset info:
runfilz = dir(sprintf('%s%cHierarchicalPriors_subject%02d_block*',indirPTB, filesep, isub));
dataPTB = [];
for irun = 1:length(runfilz)

    load(fullfile(indirPTB,runfilz(irun).name),'output');

    if irun == 1
        dataPTB = output;
    else
        dataPTB = [dataPTB output];
    end
end

% clean up data: remove data after test onset for catch trials, and remove bad segments
% remove data after test onset for each catch trial
catchtrials = find(dataMEG.trialinfo>100);
catch_time = extractfield(dataPTB(catchtrials),'time_catch');
catch_time = dsearchn(dataMEG.time',catch_time');% change from time to samples
for icatch = 1:length(catch_time)
    dataMEG.trial(catchtrials(icatch),:,catch_time(icatch):end) = NaN;
end

badfile = sprintf('%s%cbadtrials_%s.mat',indirMEG,filesep,subjectKEY{isub,2});
load(badfile,'allBadTrials');

% compare triggers as check
if sum(extractfield(dataPTB,'trigger') ~= dataMEG.trialinfo') > 0
    error('trigger values not matching between MEG and PTB');
end

allBadTrials = find(allBadTrials);
% remove bad segments
for ibadtrial = 1:length(allBadTrials)

    badID = contains(extractfield(BSTevents{allBadTrials(ibadtrial)},'label'),'bad');
    badID = find(badID);

    for ibadtype = 1:length(badID)

        bad_times = BSTevents{allBadTrials(ibadtrial)}(badID(ibadtype)).times;

        for ibadseg = 1:size(bad_times,2)

            startID = dsearchn(dataMEG.time',bad_times(1,ibadseg));
            endID = dsearchn(dataMEG.time',bad_times(2,ibadseg));

            dataMEG.trial(allBadTrials(ibadtrial),:,startID:endID) = NaN;
        end

    end

end

% now change catch trial triggers to normal values
dataMEG.trialinfo(catchtrials) = dataMEG.trialinfo(catchtrials)-100;

% for some dissimilarity measures (like correlation as used here) we can now average over trials, which highly reduces storage
if dissimilarity == 1
    triggers = unique(dataMEG.trialinfo);
    temp = zeros(triggers(end),size(dataMEG.trial,2),size(dataMEG.trial,3));
    for itrig = triggers'
        temp(itrig,:,:) = squeeze(mean(dataMEG.trial(dataMEG.trialinfo == itrig,:,:),'omitnan'));
    end
    dataMEG.trial = temp(triggers,:,:);
    dataMEG.trialinfo = triggers;
    dataMEG.sampleinfo = [];
    clear temp
end

% MNE source reconstruction using sensor data and inversion kernel created in Brainstorm
sourcefile = dir(fullfile(indirSource,'*MEG_GRAD_MEG_MAG_KERNEL_230823*'));
fn2load = sprintf('%s%c%s',indirSource,filesep,sourcefile.name);
kernel = load(fn2load);
chanIdx = kernel.GoodChannel;
cortex = kernel.SurfaceFile(6:end-4);
kernel = kernel.ImagingKernel;

atlasfile = dir(fullfile(indirAtlas,'*cortex*'));
atlasID = contains(extractfield(atlasfile,'name'),cortex);
fn2load = sprintf('%s%c%s',indirAtlas,filesep,atlasfile(atlasID).name);
atlas = load(fn2load);
atlasIdx = contains(extractfield(atlas.Atlas,'Name'),atlas2use);
atlas = atlas.Atlas(atlasIdx).Scouts;

% do some atlas fixing:
% 1) in Schaefer there is also a background
% 2) in HCP some of the L and R to indicate hemisphere behind the parcel names are gone
idx2remove = contains(extractfield(atlas,'Label'),'Background');
atlas(idx2remove) = [];

for iparcel = 1:length(atlas)

    if ~contains(atlas(iparcel).Label,' L') && ~contains(atlas(iparcel).Label,' R')

        if strcmp(atlas(iparcel).Region(1),'L')
            atlas(iparcel).Label = [atlas(iparcel).Label ' L'];
        elseif strcmp(atlas(iparcel).Region(1),'R')
            atlas(iparcel).Label = [atlas(iparcel).Label ' R'];
        end

    end

end

% select only channels that MNE was applied to
dataMEG.label = dataMEG.label(chanIdx);
dataMEG.trial = dataMEG.trial(:,chanIdx,:);

% select ROI
ROIsel = [];
ROIsel.names = ROIdefinition.names(iROI);
ROIsel.parcels = ROIdefinition.parcels(iROI);
new = ROIsel;

% divide in left and right hemisphere
for isource = 1:length(ROIsel.names)

    for iparcel = 1:length(ROIsel.parcels{isource})

        new.parcels{isource}{iparcel+length(ROIsel.parcels{isource})} = [ROIsel.parcels{isource}{iparcel} ' R'];
        new.parcels{isource}{iparcel} = [ROIsel.parcels{isource}{iparcel} ' L'];

    end

end
ROIsel = new; clear new;

% select parcels for current ROI, select kernel for those, and apply to data
vertices = [];
for iparcel = 1:length(ROIsel.parcels{isource})
    parcelidx = find(strcmp(extractfield(atlas,'Label'),ROIsel.parcels{isource}{iparcel}),1);

    if isempty(parcelidx)
        error(['cannot find matching label in atlas for ' ROIsel.parcels{isource}{iparcel}]);
    end

    vertices = [vertices atlas(parcelidx).Vertices(:).'];

end

ROIsel.vertices{isource} = vertices;

% select inversion kernel for this ROI
kernel = kernel(vertices,:);

% multiple sensor level data with inversion kernel to get source level data
sourcedata = zeros(size(dataMEG.trial,1),size(kernel,1),size(dataMEG.trial,3));
for itrial = 1:size(dataMEG.trial,1)
    temp = squeeze(dataMEG.trial(itrial,:,:));
    sourcedata(itrial,:,:) = kernel*temp;
end
dataMEG.trial = single(sourcedata);clear sourcedata temp kernel
dataMEG.label = vertices';

%% load model RDMs
load([indirMOD filesep 'HierarchicalPriors_dynRDM_graysmooth'],'RDMgraysmooth');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_OFmag'],'RDMoptflow_mag');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_OFdir'],'RDMoptflow_dir');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_posturedep'],'RDMposture_dep');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_postureinvar'],'RDMposture_invar');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_motiondep'],'RDMmotion_dep');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_motioninvar'],'RDMmotion_invar');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_accdep'],'RDMacc_dep');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_accinvar'],'RDMacc_invar');

dynRDM{1} = single(RDMgraysmooth);
dynRDM{2} = single(RDMoptflow_mag);
dynRDM{3} = single(RDMoptflow_dir);
dynRDM{4} = single(RDMposture_dep);
dynRDM{5} = single(RDMposture_invar);
dynRDM{6} = single(RDMmotion_dep);
dynRDM{7} = single(RDMmotion_invar);
dynRDM{8} = single(RDMacc_dep);
dynRDM{9} = single(RDMacc_invar);

% clear up workspace
clear RDMgraysmooth RDMoptflow_mag RDMoptflow_dir RDMposture_dep RDMposture_invar RDMmotion_dep RDMmotion_invar RDMacc_dep RDMacc_invar

% downsample models if needed
if fsNew ~= 100
    for imodel = 1:length(dynRDM)
        decimatefactor = (size(dynRDM{imodel},3)/5)/fsNew;
        dynRDM{imodel} = dynRDM{imodel}(:,:,1:decimatefactor:end,1:decimatefactor:end);
    end
end

% time definitions
fsKin = fsNew;% sample rate kinematic recordings
tKin = 0:1/fsKin:5-1/fsKin;
tNeural = dataMEG.time;
oldtimelen = length(tKin)/fsKin;
latencytime = -maxlatency:1/fsNew:maxlatency;% for latency plots
lagIDs = - maxlatency*fsNew:maxlatency*fsNew;% for averaging over diagonals in the 2D dRSA matrix

% loop over conditions
for icon = conditions

    % I forgot to store the proper scramble indices for the first 3 subjects so we can't use the scramble condition, but the other 2 conditions are usable
    if isub < 4 && icon == 3
        exit
    end

    fn2save = sprintf('%s%cSUB%02d_%s_%dHz_%s_dissimilarity%d_smMEG%d_smRDMneu%d_smRDMmod%d',...
        outdir, filesep, isub, conditionnames{icon}, fsNew, ROIdefinition.names{iROI}{:}, dissimilarity, smoothNeuralData, smoothNeuralRDM, smoothModelRDM);

    if exist([fn2save '.mat'],'file')
        continue
    end

    % Select trials for current condition
    dataMEGcon = dataMEG;
    current_triggers = 1:14;
    current_triggers = current_triggers+20*(icon-1);
    trials2select = ismember(dataMEG.trialinfo,current_triggers);
    dataMEGcon.trial = dataMEGcon.trial(trials2select,:,:);
    dataMEGcon.trialinfo = current_triggers;

    %% load and prepare eyetracker RDMs
    fn= sprintf('%s%cSUB%02d_dynRDM_eyeTracker',indirEYE, filesep, isub);
    load(fn,['RDMeye' upper(conditionnames{icon})]);
    dynRDM{10} = eval(['single(RDMeye' upper(conditionnames{icon}) ')']);
    clear(['RDMeye' upper(conditionnames{icon})]);

    % downsample models if needed, because models are at 100 Hz
    decimatefactor = (size(dynRDM{10},3)/oldtimelen)/fsNew;
    dynRDM{10} = dynRDM{10}(:,:,1:decimatefactor:end,1:decimatefactor:end);

    %% piece-wise scramble model time lines according to individual subject scramble indices (i.e., as they saw the video)
    if icon == 3
        % load individual subject scrambled video indices
        fn = sprintf('%s%cHierarchicalPriors_subject%02d_scramindices',indirPTB,filesep,isub);
        load(fn,'samples_all','randID_all');

        % loaded indices are at 120 Hz at which videos were displayed, so convert to cfg.fsNew so it matches the model RDMs
        % and create vector of indices to scramble the model RDMs with
        scrambledIDs = zeros(nstim,length(tKin));
        for istim = 1:nstim
            temp = round(samples_all{istim} * (fsNew/120));
            temp(temp==0) = 1;% if fsNew < 120/2, the first index will be zero, which is erroneous for an index

            % sometimes end of segment n overlaps with beginning of segment n+1, so fix that
            segment2fix = find(temp(2:end,1) ~= temp(1:end-1,2)+1);
            temp(segment2fix,2) = temp(segment2fix+1,1)-1;

            %
            temp2 = [];
            for isegment = 1:size(temp,1)
                start_end = temp(randID_all{istim}(isegment),:);
                temp2 = [temp2 start_end(1):start_end(2)];% add all new indices in single vector
            end

            % sanity check
            if any(sort(temp2)~=1:length(temp2))
                error('error in scramble indices');
            end

            scrambledIDs(istim,:) = temp2;% store new indices

        end
        clear temp temp2

        % now scramble the model RDMs
        for imodel = 1:length(dynRDM)-1% don't scramble the last RDM which is based on eyetracker data
            for istim1 = 1:nstim
                for istim2 = 1:nstim

                    dynRDM{imodel}(istim1,istim2,:,:) = dynRDM{imodel}(istim1,istim2,scrambledIDs(istim1,:),scrambledIDs(istim2,:));

                end
            end
        end

    end

    % ----------------------------------------------------------------------------------------------------
    % ----------------------------------------------------------------------------------------------------
    % create matrix with times for randomly shuffling sequence onsets over X iterations
    % shuffleTime = start times size [iterations X stimuli];
    rng((round(sum(clock))+icon*10000+isub*1000+iROI)*10000);% make sure rand numbers are different for each subject, each ROI, each condition, and each time this script is ran at a different time and date
    onsetIDs = 0.2+round((rand(iternum,nstim)*(oldtimelen-stimlen-0.2))./(1/fsNew)).*(1/fsNew);% make sure it starts from 200 msec (no earlier prediction possible), and make sure it's a multiple of 1/fsNew, such that exact same time point is picked for neural and kinematic data

    dRSAitercumsum = zeros(length(latencytime),length(models2test),'single');
    for iter = 1:iternum

        % determine indices for stimulus-specific re-alignment for current iteration
        onsetIDmodel = dsearchn(tKin',onsetIDs(iter,:)')';
        onsetIDneural = dsearchn(tNeural',onsetIDs(iter,:)')';

        %% prepare neural RDM for current iteration
        dataMEGiter = dataMEGcon;

        % re-align neural data to stimulus-specific temporal jitter
        temptrial = zeros(nstim,length(dataMEGiter.label),fsNew*stimlen,'single');

        for istim = 1:nstim

            tNewID = onsetIDneural(istim):onsetIDneural(istim)+stimlen*fsNew-1;

            temptrial(istim,:,:) = squeeze(dataMEGiter.trial(istim,:,tNewID));

        end
        dataMEGiter.trial = temptrial; clear temptrial
        dataMEGiter.time = 0:1/fsNew:stimlen-1/fsNew;
        dataMEGiter.trialinfo = 1:nstim;

        % centre data around conditions (recommended by some incl. e.g., cosmomvpa)
        dataMEGiter.trial = dataMEGiter.trial - repmat(mean(dataMEGiter.trial,'omitnan'),size(dataMEGiter.trial,1),1,1);

        % compute neural RDM using correlation, this short section could be replaced by other dissimilarity measure
        MATxTIME = zeros(nstim,nstim,length(dataMEGiter.time),'single');
        for itime = 1:length(dataMEGiter.time)
            MATxTIME(:,:,itime) = 1-corr(squeeze(dataMEGiter.trial(:,:,itime))');
        end% time loop

        % smooth neural RDM across time
        if smoothNeuralRDM > 0
            for i=1:size(MATxTIME,1)
                MATxTIME(i,:,:) = ft_preproc_smooth(squeeze(MATxTIME(i,:,:)),smoothNeuralRDM);
            end
        end

        % vectorize the matrices
        neuralRDM = zeros((nstim*nstim-nstim)/2,size(MATxTIME,3),'single');
        for itime = 1:size(MATxTIME,3)
            MAT = squeeze(MATxTIME(:,:,itime));
            neuralRDM(:,itime) = squareform(MAT);
        end

        %% prepare model RDM for current iteration
        % subsample and re-align model RDMs for current iteration 
        % (in practice you're extracting a smaller [1 stim x 1 stim x 3 sec x 3 sec] square from the larger [14 x 14 x 5 sec x 5 sec] dynamic RDMs),
        % and then take the diagonal of that smaller square. 
        modelRDMsquare = zeros(length(dynRDM),nstim,nstim,fsNew*stimlen,'single');
        for istim1 = 1:nstim
            for istim2 = 1:nstim

                tNewID1 = onsetIDmodel(istim1):onsetIDmodel(istim1)+stimlen*fsNew-1;
                tNewID2 = onsetIDmodel(istim2):onsetIDmodel(istim2)+stimlen*fsNew-1;

                for imodel = 1:length(dynRDM)

                    temptrial = squeeze(dynRDM{imodel}(istim1,istim2,tNewID1,tNewID2));

                    temptrial = diag(temptrial);

                    modelRDMsquare(imodel,istim1,istim2,:) = temptrial;

                end% model loop
            end% second stimulus loop
        end% first stimulus loop

        % for view-invariant models take average of (stim1*stim2trans) and (stim1trans*stim2), which is not exactly the same because procrustes is not symmetrical
        for imodel = [5 7 9]
            for itime = 1:size(modelRDMsquare,4)
                temp = squeeze(modelRDMsquare(imodel,:,:,itime));
                modelRDMsquare(imodel,:,:,itime) = (temp + temp')/2;
            end% time loop
        end% model loop

        %% some more model RDM preparation
        % extract vector RDM from triangle of square RDM
        modelRDM = zeros(size(modelRDMsquare,1),size(modelRDMsquare,4),(nstim*nstim-nstim)/2,'single');
        for imodel = 1:size(modelRDMsquare,1)
            for itime = 1:size(modelRDMsquare,4)
                modelRDM(imodel,itime,:) = squareform(tril(squeeze(modelRDMsquare(imodel,:,:,itime)),-1));
            end
        end

        % smooth model RDM across time
        if smoothModelRDM% && icon ~=3% if scrambled, smoothing already happened above pre-scrambling, but because that's slower, better do it here for other conditions
            for i=1:size(modelRDM,3)
                modelRDM(:,:,i) = ft_preproc_smooth(squeeze(modelRDM(:,:,i)),smoothModelRDM);
            end
        end

        % centre neuralRDM per time point
        neuralRDM = neuralRDM - repmat(mean(neuralRDM,'omitnan'),size(neuralRDM,1),1);

        % standardize neuralRDM across all time points at once to keep temporal structure
        neuralRDM = neuralRDM ./ std(neuralRDM(:));

        % centre and standardize modelRDM
        for imodel = 1:size(modelRDM,1)
            temp = squeeze(modelRDM(imodel,:,:))';

            % centre modelRDM per model and per time point
            temp = temp - repmat(mean(temp,'omitnan'),size(temp,1),1);

            % standardize modelRDM across all time points at once to keep temporal structure
            modelRDM(imodel,:,:) = temp' ./ std(temp(:));

        end

        %% dynamic RSA
        if similarity == 0% simple correlation
            dRSA = zeros(size(neuralRDM,2),size(neuralRDM,2),length(models2test),'single');
            for imodel = 1:length(models2test)% only test for our standard 7

                % neural - model correlation
                dRSA(:,:,imodel) = corr(squeeze(modelRDM(imodel,:,:))',neuralRDM);

            end

        elseif similarity == 1% use regression to regress out principal components of other models

            % which models to regress out:
            models2regressout = [2:10; 1 3:10; 1 2 4:10; 1:3 5:10; 1:4 6:10; 1:5 7:10; 1:6 8:10; 1:7 9:10; 1:8 10; 1:9];
            models2regressout = models2regressout(models2test,:);

            load(fullfile(outdir,'..',['regressionBorder_smRDM30msec_' num2str(fsNew) 'Hz_' conditionnames{icon}]),'regborder');
            regborder.subinvarmods(10) = regborder.subvarmods(isub == subjects);
            regborder = regborder.subinvarmods(models2test);

            dRSA = zeros(size(neuralRDM,2),size(neuralRDM,2),length(models2test),'single');
            for imodel = 1:length(models2test)% models to test

                for itime = 1:size(modelRDM,2)% model time

                    % regressor selection: 1. only regress out in the -1.5 to 1.5 lag interval, 2. discard ones outside of video time
                    regidx = itime-fsNew*1.5:itime+fsNew*1.5;
                    regidx(logical((regidx<1) + (regidx>size(modelRDM,2)))) = [];

                    % selection of indices at which to regress out the model itself, based on regression borders, but don't bother if distance is larger than 1 sec
                    if regborder(imodel) < fsNew
                        test2regidx = [itime-fsNew*1.5:itime-regborder(imodel) itime+regborder(imodel):itime+fsNew*1.5];
                        test2regidx(logical((test2regidx<1) + (test2regidx>size(modelRDM,2)))) = [];
                    else
                        test2regidx = [];
                    end

                    % select the model RDM at the time point to test on this iteration
                    Xtest = squeeze(modelRDM(models2test(imodel),itime,:));

                    % select the to-be-tested model RDMs but at distant time points from the to-be-tested time point (with minimum distance based on regression border)
                    Xtest2regressout = squeeze(modelRDM(models2test(imodel),test2regidx,:))';

                    % then use PCA to reduce dimensionality by including only components that explain at least 0.1% of total variance
                    [~, score, ~, ~, exp, ~] = pca(Xtest2regressout);
                    imax = sum(exp>.1);
                    Xtest2regressout = score(:,1:imax);

                    % same as above, but now for regressing out all other models across the entire latency window defined by regidx
                    Xregressout = zeros(size(modelRDM,3),500);% add some zeros for initialization
                    for ireg = 1:nnz(models2regressout(imodel,:))

                        temp = squeeze(modelRDM(models2regressout(imodel,ireg),regidx,:))';
                        [~, score, ~, ~, exp, ~] = pca(temp);
                        imax = sum(exp>.1);% only components with minimum variance of X%
                        score = score(:,1:imax);
                        Xregressout(:,nnz(Xregressout(1,:))+1:nnz(Xregressout(1,:))+size(score,2)) = score;

                    end
                    Xregressout = Xregressout(:,1:nnz(Xregressout(1,:)));% remove zeros
                    Xregressout = [Xregressout Xtest2regressout];% combine all RDMs to-be-regressed out from to-be-tested model itself and all other models

                    % Xregressout needs to still be standardized again
                    Xregressout = Xregressout/std(Xregressout(:));

                    % prepare regressors and response variable for PCR
                    X = [Xtest Xregressout];% regressors, with our interest in the beta-weight of Xtest
                    Y = neuralRDM;% response variable

                    % run principal component regression (PCR)
                    [PCALoadings,PCAScores] = pca(X,'NumComponents',nPCRcomps);% compute PCA components of original regressors
                    betaPCR = PCAScores\(Y);% regression on PCA components instead of original regressors
                    betaPCR = PCALoadings*betaPCR;% back-project beta-weights of PCA components to space of original regressors by using the PCA loadings

                    % extract only weights for Xtest and ignore weights for Xregressout
                    dRSA(itime,:,imodel) = betaPCR(1,:);

                end% ibin1 loop

            end% model loop

        end% similarity if statement

        % average over diagonals in dRSA matrix and combine data from all iterations
        for imodel = 1:size(dRSA,3)

            temp = spdiags(squeeze(dRSA(:,:,imodel)), lagIDs);
            temp(temp == 0) = nan;
            dRSAitercumsum(:,imodel) = dRSAitercumsum(:,imodel) + mean(temp,'omitnan')';

        end% model loop
        clear dRSA temp

        % store a progress file so we can see the progress after every 100 iterations even if the script runs on the server
        if mod(iter,100) == 0
            fnprog = sprintf('%s%cPROGRESS_SUB%02d_%s.mat',outdir, filesep, isub, ROIdefinition.names{iROI}{:});
            save(fnprog,'iter')
        end

    end% iteration loop

    % average over iterations
    dRSAall = dRSAitercumsum / iternum;
    clear dRSAitercumsum;

    save(fn2save,'dRSAall','ROIsel');

    % remove progress file
    delete(fnprog);

end% loop over conditions

end