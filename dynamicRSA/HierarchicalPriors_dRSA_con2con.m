function HierarchicalPriors_dRSA_con2con(parms,isub,iROI,~,~,~)

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
indirPTB = fullfile(rootdir,'data','MEG','PTBoutput');
indirSource = fullfile(rootdir,'data','MEG','brainstorm_database','HierarchicalPriors','data',subjectKEY{isub,2},subjectKEY{isub,2});
indirAtlas = fullfile(rootdir,'data','MEG','brainstorm_database','HierarchicalPriors','anat',subjectKEY{isub,2});

% create directory for results
simstring = {'corr_con2con' ['pcr_con2con_' num2str(nPCRcomps) 'comps'],['corr_con2con_partout' num2str(mod2part)]};
simstring = simstring{similarity+1};
outdir = sprintf('%s%sresults%ssource_%s%sRSA%s%s_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,simstring,iternum,ceil(stimlen));
if ~exist(outdir,'dir')
    mkdir(outdir);
end

load(fullfile(rootdir,'code','dRSA','ROIdefinitions'),'ROIdefinition');
conditionnames = {'normal','invert','scram'};

fn2save = sprintf('%s%cSUB%02d_%s2%s_%dHz_%s_dissimilarity%d_smMEG%d_smRDMneu%d',...
    outdir, filesep, isub, conditionnames{con2con(1)}, conditionnames{con2con(2)}, fsNew, ROIdefinition.names{iROI}{:}, dissimilarity, smoothNeuralData, smoothNeuralRDM);

if exist([fn2save '.mat'],'file')
    return
end

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
end

% For catch trials, find test onset and change all data after test onset to NaNs,
% this way, when averaging over trials using mean with the omitnan flag,
% only the pre-test interval is used
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

% for some dissimilarity measures we can now average over trials, which highly reduces storage
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

kernel = kernel(vertices,:);
sourcedata = zeros(size(dataMEG.trial,1),size(kernel,1),size(dataMEG.trial,3));
for itrial = 1:size(dataMEG.trial,1)
    temp = squeeze(dataMEG.trial(itrial,:,:));
    sourcedata(itrial,:,:) = kernel*temp;
end
dataMEG.trial = single(sourcedata);clear sourcedata temp kernel
dataMEG.label = vertices';

% time definitions
fsKin = fsNew;% sample rate kinematic recordings
tKin = 0:1/fsKin:5-1/fsKin;
tNeural = dataMEG.time;
oldtimelen = length(tKin)/fsKin;
latencytime = -maxlatency:1/fsNew:maxlatency;% for latency plots
lagIDs = - maxlatency*fsNew:maxlatency*fsNew;% for averaging over diagonals in the 2D dRSA matrix

%% load model RDMs
load([indirMOD filesep 'HierarchicalPriors_dynRDM_graysmooth'],'RDMgraysmooth');
dynRDM{1} = single(RDMgraysmooth);
if mod2part == 2
    load([indirMOD filesep 'HierarchicalPriors_dynRDM_OFmag'],'RDMoptflow_mag');
    dynRDM{2} = single(RDMoptflow_mag);
end

% clear up workspace
clear RDMgraysmooth RDMoptflow_mag

% downsample or upsample models if needed
if fsNew < 100
    for imodel = 1:length(dynRDM)
        decimatefactor = (size(dynRDM{imodel},3)/5)/fsNew;
        dynRDM{imodel} = dynRDM{imodel}(:,:,1:decimatefactor:end,1:decimatefactor:end);
    end
elseif fsNew > 100
    for imodel = 1:length(dynRDM)
        fsOld = 100;
        tOld = 0:1/fsOld:5-1/fsOld;
        tempall = zeros(nstim,nstim,length(tKin),length(tKin));
        for istim1 = 1:nstim
            for istim2 = 1:nstim
                temp = squeeze(dynRDM{imodel}(istim1,istim2,:,:));
                % bilinear interpolation
                temp = interp1(tOld,temp,tKin,'pchip','extrap');
                temp = interp1(tOld,temp',tKin,'pchip','extrap')';
                tempall(istim1,istim2,:,:) = temp;
            end
        end
        dynRDM{imodel} = tempall;
        clear tempall temp
    end
end

% select trials for 'normal' condition as outcome variable (MEGout and RDMout), or 'neural' in standard dRSA
MEGout = dataMEG;
current_triggers = 1:14;
current_triggers = current_triggers+20*(con2con(1)-1);
trials2select = ismember(dataMEG.trialinfo,current_triggers);
MEGout.trial = dataMEG.trial(trials2select,:,:);
MEGout.trialinfo = current_triggers;

% select trials for 'invert' condition as regressor variable (MEGreg and RDMreg), or 'model' in standard dRSA
MEGreg = dataMEG;
current_triggers = 1:14;
current_triggers = current_triggers+20*(con2con(2)-1);
trials2select = ismember(dataMEG.trialinfo,current_triggers);
MEGreg.trial = dataMEG.trial(trials2select,:,:);
MEGreg.trialinfo = current_triggers;

% ----------------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------
% create matrix with times for randomly shuffling sequence onsets over X iterations
% shuffleTime = start times size [iterations X stimuli];
rng((round(sum(clock))+isub*1000+iROI)*10000);% make sure rand numbers are different for each subject, each ROI, each condition, and each time this script is ran at a different time and date
onsetIDs = 0.2+round((rand(iternum,nstim)*(oldtimelen-stimlen-0.2))./(1/fsNew)).*(1/fsNew);% make sure it starts from 200 msec (no earlier prediction possible), and make sure it's a multiple of 1/fsNew, such that exact same time point is picked for neural and kinematic data

dRSAitercumsum = zeros(length(latencytime),1,'single');
for iter = 1:iternum

    % determine indices for stimulus-specific re-alignment for current iteration
    onsetIDneural = dsearchn(tNeural',onsetIDs(iter,:)')';
    onsetIDmodel = dsearchn(tKin',onsetIDs(iter,:)')';

    %% prepare neural RDM for current iteration
    for icon2con = 1:2

        if icon2con == 1
            dataMEGiter = MEGout;
        else
            dataMEGiter = MEGreg;
        end

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

        % compute neural RDM using correlation
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
        neuralRDM = zeros((nstim*nstim-nstim)/2,length(dataMEGiter.time),'single');
        for itime = 1:size(MATxTIME,3)
            MAT = squeeze(MATxTIME(:,:,itime));
            neuralRDM(:,itime) = squareform(MAT);
        end

        % centre neuralRDM per time point
        neuralRDM = neuralRDM - repmat(mean(neuralRDM,'omitnan'),size(neuralRDM,1),1);

        % standardize neuralRDM across all time points at once to keep temporal structure
        neuralRDM = neuralRDM ./ std(neuralRDM(:));

        if icon2con == 1
            RDMout = neuralRDM;
        else
            RDMreg = neuralRDM;
        end

    end% con2con loop

    %% prepare model RDM for current iteration
    % subsample and re-align model RDMs for current iteration
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

    % centre and standardize
    for imodel = 1:size(modelRDM,1)
        temp = squeeze(modelRDM(imodel,:,:))';

        % centre modelRDM per model and per time point
        temp = temp - repmat(mean(temp,'omitnan'),size(temp,1),1);

        % standardize modelRDM across all time points at once to keep temporal structure
        modelRDM(imodel,:,:) = temp' ./ std(temp(:));
    end

    %% dynamic RSA
    if similarity == 0% neural - neural correlation

        dRSA = corr(RDMreg,RDMout);

    elseif similarity == 2% partial correlation, controlling for low-level visual model(s)

        if size(modelRDM,1) == 1% if only 1 model gets partialled out
            modelRDM = squeeze(modelRDM)';
        else% concatenate if multiple get partialled out
            modelRDM = reshape(modelRDM,size(modelRDM,1)*size(modelRDM,2),size(modelRDM,3))';
        end

        % first reduce dimensionality of model RDM
        [~, score, ~, ~, exp, ~] = pca(modelRDM,'NumComponents',75);% should be well below number of observations (i.e., 91)
        imax = find(cumsum(exp)>95,1);% take first x components that together explain at least 95% variance
        modelRDM = score(:,1:imax);

        dRSA = partialcorr(RDMreg,RDMout,modelRDM);

    end% similarity if loop

    % average over diagonals in dRSA matrix and combine data from all iterations
    temp = spdiags(dRSA, lagIDs);
    temp(temp == 0) = nan;
    dRSAitercumsum = dRSAitercumsum + mean(temp,'omitnan')';

    clear dRSA temp

end% iteration loop

% average over iterations
dRSAall = dRSAitercumsum / iternum;
clear dRSAitercumsum;

save(fn2save,'dRSAall','ROIsel');

end