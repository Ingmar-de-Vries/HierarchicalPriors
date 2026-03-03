function HierarchicalPriors_dRSA_simulations(parms,isub,~,~,~,isim)

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
indirMOD = fullfile(rootdir,'data','modelRDMs');
indirEYE = fullfile(rootdir,'data','MEG','ELpreprocessed');
indirPTB = fullfile(rootdir,'data','MEG','PTBoutput');

% create directory for results
simstring = {'corr' ['pcaANDpcr_' num2str(nPCRcomps) 'comps']};
simstring = simstring{similarity+1};
outdir = sprintf('%s%sresults%ssource_%s%sRSA%ssimulations%s%s_%dlag_%diterations_%dsec',...
    rootdir,filesep,filesep,atlas2use,filesep,filesep,filesep,simstring,lag*1000,iternum*length(subjects),ceil(stimlen));
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% time definitions
fsKin = fsNew;% sample rate kinematic recordings
tKin = 0:1/fsKin:5-1/fsKin;
oldtimelen = length(tKin)/fsKin;
latencytime = -maxlatency:1/fsNew:maxlatency;% for latency plots
lagIDs = - maxlatency*fsNew:maxlatency*fsNew;% for averaging over diagonals in the 2D dRSA matrix

%% set weights for model to implant
if lag ~= 0
    error('IF YOU WANT TO INTRODUCE A LAG IN THE SIMULATED MODEL, YOU NEED TO USE THE PREVIOUS INTERPOLATION METHOD THAT IS COMMENTED NOW');
end

simWeights = zeros(1,length(latencytime));
simWeights(dsearchn(latencytime',lag)) = 1;
simWeights = repmat(simWeights,99,1);
if isim == 99
    simWeights = zeros(size(simWeights));
end
zerolag = ceil(size(simWeights,2)./2);% for the simulations

%% load model RDMs
load([indirMOD filesep 'HierarchicalPriors_dynRDM_graysmooth'],'RDMgraysmooth');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_OFdir'],'RDMoptflow_dir');
load([indirMOD filesep 'HierarchicalPriors_dynRDM_OFmag'],'RDMoptflow_mag');
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

conditionnames = {'normal','invert','scram'};
for icon = conditions

    % I forgot to store the proper scramble indices for the first 3 subjects so we can't use the scramble condition, but the other 2 conditions are usable
    if isub < 4 && icon == 3
        exit
    end

    fn2save = sprintf('%s%cSUB%02d_%s_implant%d_%dHz_dissimilarity%d_smMEG%d_smRDMneu%d_smRDMmod%d',...
        outdir, filesep, isub, conditionnames{icon}, isim, fsNew, dissimilarity, smoothNeuralData, smoothNeuralRDM, smoothModelRDM);

    %% load and prepare eyetracker RDMs
    fn= sprintf('%s%cSUB%02d_dynRDM_eyeTracker',indirEYE, filesep, isub);
    load(fn,['RDMeye' upper(conditionnames{icon})]);
    dynRDM{10} = eval(['single(RDMeye' upper(conditionnames{icon}) ')']);
    clear(['RDMeye' upper(conditionnames{icon})]);

    % downsample models if needed
    if fsNew ~= 100
        for imodel = 1:length(dynRDM)
            decimatefactor = (size(dynRDM{imodel},3)/5)/fsNew;
            dynRDM{imodel} = dynRDM{imodel}(:,:,1:decimatefactor:end,1:decimatefactor:end);
        end
    end

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
    rng((round(sum(clock))+icon*10000+isub*1000+isim)*10000);% make sure rand numbers are different for each subject, each ROI, each condition, and each time this script is ran at a different time and date
    onsetIDs = 0.2+round((rand(iternum,nstim)*(oldtimelen-stimlen-0.2))./(1/fsNew)).*(1/fsNew);% make sure it starts from 200 msec (no earlier prediction possible), and make sure it's a multiple of 1/fsNew, such that exact same time point is picked for neural and kinematic data

    % subsampling loop
    dRSAitercumsum = zeros(length(latencytime),length(modelnames),'single');
    for iter = 1:iternum

        % determine indices for stimulus-specific re-alignment for current iteration
        onsetIDmodel = dsearchn(tKin',onsetIDs(iter,:)')';

        %% prepare and implant model RDM into randomized neural RDM to create simulated RDM
        % randomize neural RDM as basis for simulated neural RDM, random uniformly distributed numbers between 0 and 2
        neuralRDM = rand((nstim*nstim-nstim)./2,fsNew*stimlen);

        % smooth random neural RDM across time
        if smoothNeuralRDM
            neuralRDM = ft_preproc_smooth(neuralRDM,smoothNeuralRDM);
        end

        % create model RDM that will be implanted into the random neural RDM according to [lag X model] weight matrix in simWeights
        RDM2implant = zeros(nstim,nstim,fsNew*stimlen,'single');
        for istim1 = 1:nstim
            for istim2 = 1:nstim

                % NOTE: IF YOU WANT TO INTRODUCE A LAG IN THE SIMULATED MODEL, YOU NEED TO USE THE PREVIOUS INTERPOLATION METHOD THAT IS COMMENTED BELOW
                % determine new time which automatically realigns resampled data
                % tKin4extrap1 = [-tstart(istim1)/cfg.newsampfreq:1/cfg.newsampfreq:tKin(1) tKin(2:end-1) tKin(end):1/cfg.newsampfreq:tKin(end)+(cfg.subsample(2)*cfg.newsampfreq-tstart(istim1))/cfg.newsampfreq];
                % tKin4extrap2 = [-tstart(istim2)/cfg.newsampfreq:1/cfg.newsampfreq:tKin(1) tKin(2:end-1) tKin(end):1/cfg.newsampfreq:tKin(end)+(cfg.subsample(2)*cfg.newsampfreq-tstart(istim2))/cfg.newsampfreq];

                tNewID1 = onsetIDmodel(istim1):onsetIDmodel(istim1)+fsNew*stimlen-1;
                tNewID2 = onsetIDmodel(istim2):onsetIDmodel(istim2)+fsNew*stimlen-1;

                lags = find(simWeights(isim,:) ~= 0);
                for ilag = 1:length(lags)% loop over different lags at which current model is implanted

                    currentlag = lags(ilag);
                    weight = simWeights(isim,currentlag);
                    lag2implant = currentlag - zerolag;
                    lag2implant = lag2implant / fsNew;%change from samples to seconds

                    %                 temptrial = interp1(tKin,squeeze(dynRDM{implant}(istim1,istim2,:,:)),tKin4extrap1-lag2implant,'pchip','extrap');
                    %                 temptrial = interp1(tKin,temptrial',tKin4extrap2-lag2implant,'pchip','extrap')';

                    temptrial = squeeze(dynRDM{isim}(istim1,istim2,tNewID1+lag2implant,tNewID2+lag2implant));

                    % keep extrapolated values within boundaries:
                    %                 maxval = max(max(temptrial(tstart(istim1):tstart(istim1)+length(modelID),tstart(istim2):tstart(istim2)+length(modelID))));
                    %                 minval = min(min(temptrial(tstart(istim1):tstart(istim1)+length(modelID),tstart(istim2):tstart(istim2)+length(modelID))));
                    %                 temptrial(temptrial > maxval) = maxval;
                    %                 temptrial(temptrial < minval) = minval;

                    temptrial = diag(temptrial);
                    RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + weight*temptrial;
                end% lag loop

            end% second stimulus loop
        end% first stimulus loop

        % for view-invariant models take average of (stim1*stim2trans) and (stim1trans*stim2), which is not exactly the same because procrustes is not symmetrical
        if isim == 5 || isim == 7 || isim == 9
            for itime = 1:size(RDM2implant,3)
                temp = squeeze(RDM2implant(:,:,itime));
                RDM2implant(:,:,itime) = (temp + temp')/2;
            end% time loop
        end

        % extract vector RDM from triangle of square RDM
        temp = zeros((nstim*nstim-nstim)/2,size(RDM2implant,3));
        for itime = 1:size(RDM2implant,3)
            temp(:,itime) = squareform(tril(squeeze(RDM2implant(:,:,itime)),-1));
        end
        RDM2implant = temp;clear temp

        % smooth implanted RDM across time
        if smoothNeuralRDM
            RDM2implant = ft_preproc_smooth(RDM2implant,smoothNeuralRDM);
        end

        % Implant model RDM into randomized neural RDM
        if isim ~= 99% if 99, don't implant anything
            neuralRDM = randomNeuralRDMweight*neuralRDM + RDM2implant;
        end

        % centre neuralRDM per time point
        neuralRDM = neuralRDM - repmat(mean(neuralRDM,'omitnan'),size(neuralRDM,1),1);

        % standardize neuralRDM across all time points at once to keep temporal structure
        neuralRDM = neuralRDM ./ std(neuralRDM(:));

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
        if smoothModelRDM% if scrambled, smoothing already happened above pre-scrambling, but because that's slower, better do it here for other conditions
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
        if similarity == 0% simple correlation
            dRSA = zeros(size(neuralRDM,2),size(neuralRDM,2),length(modelnames),'single');
            for imodel = 1:length(modelnames)% only test for our standard 7

                % neural - model correlation
                dRSA(:,:,imodel) = corr(squeeze(modelRDM(imodel,:,:))',neuralRDM);

            end

        elseif similarity == 1% use regression to regress out principal components of other models
            
            % which models to regress out:
            models2regressout = [2:10; 1 3:10; 1 2 4:10; 1:3 5:10; 1:4 6:10; 1:5 7:10; 1:6 8:10; 1:7 9:10; 1:8 10; 1:9];
            
            load(fullfile(outdir,'..','..',['regressionBorder_smRDM30msec_' num2str(fsNew) 'Hz_' conditionnames{icon}]),'regborder');
            regborder.subinvarmods(10) = regborder.subvarmods(isub == subjects);
            regborder = regborder.subinvarmods;
            
            dRSA = zeros(size(neuralRDM,2),size(neuralRDM,2),length(modelnames),'single');
            for imodel = 1:length(modelnames)% models to test
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
                    Xtest = squeeze(modelRDM(imodel,itime,:));

                    % select the to-be-tested model RDMs but at distant time points from the to-be-tested time point (with minimum distance based on regression border)
                    Xtest2regressout = squeeze(modelRDM(imodel,test2regidx,:))';

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

        end

        % average over diagonals in dRSA matrix and combine data from all iterations
        for imodel = 1:size(dRSA,3)

            temp = spdiags(squeeze(dRSA(:,:,imodel)), lagIDs);
            temp(temp == 0) = nan;
            dRSAitercumsum(:,imodel) = dRSAitercumsum(:,imodel) + mean(temp,'omitnan')';

        end% model loop
        clear dRSA temp

    end

    % average over iterations
    dRSAall = dRSAitercumsum / iternum;
    clear dRSAitercumsum;

    save(fn2save,'dRSAall');

end% loop over conditions

end