function HierarchicalPriors_eyetracking_asc2ft(parms,isub)

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
indir = fullfile(rootdir,'data','MEG','ELoutput');
indirTrials = fullfile(rootdir,'data','MEG','PTBoutput');

outdir = fullfile(rootdir,'data','MEG','ELpreprocessed');
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% if on cluster, only send a single subject and single ROI combination as
% job to cluster, as indicated by isub and iroi in cluster_shell.m script
if runcluster
    subjects = isub;
end

subfilz = dir(fullfile(indir,'*1.asc*'));
trialfilz = dir(fullfile(indirTrials,'*block1*'));

missedSampleNum = cell(length(subjects),8);
missedTrialNum = zeros(length(subjects),8);
for iSub = subjects
    
    runfilz = dir(fullfile(indir,['*' subfilz(iSub).name(1:end-5) '*.asc']));
    trialrunfilz = dir(fullfile(indirTrials,['*' trialfilz(iSub).name(1:end-5) '*']));
    
    if iSub == 33
        trialrunfilz(2) = [];
    end
    
    if length(runfilz) ~= length(trialrunfilz)
        error('not same amount of runs!');
    end
    
    ds = cell(length(runfilz),1);
    for iRun = 1:length(runfilz)
        
        %% get trial conditions
        trialName = fullfile(indirTrials, trialrunfilz(iRun).name);
        load(trialName,'output');
        
        %% load EL data
        fprintf('subject %d, run %d  \n',iSub,iRun);
        ascName = fullfile(indir, runfilz(iRun).name);
        
        %% load asc events to get time points
        fid = fopen(ascName, 'rt');
        aline = fread(fid, inf, 'char=>char');          % returns a single long string
        fclose(fid);
        
        aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
        aline = tokenize(aline, uint8(newline));        % split on newline
        
        % first go through aline and remove fixations, saccades, blinks, and catch trial onsets, coz we only need actual position data, and video onsets
        lines2remove = zeros(numel(aline),1);
        for i=1:numel(aline)
            tline = aline{i};
            % We don't need the automatic fixation, saccade or blink registrations, only the x,y positions and pupil dilation for each sample. So skip over these:
            if  ~isempty(regexp(tline, '^SFIX' ,'once')) || ~isempty(regexp(tline, '^EFIX' ,'once')) || ~isempty(regexp(tline, '^SSACC', 'once')) || ~isempty(regexp(tline, '^ESACC', 'once')) || ~isempty(regexp(tline, '^SBLINK', 'once')) || ~isempty(regexp(tline, '^EBLINK', 'once')) || ~isempty(regexp(tline, 'ERROR', 'once'))
                lines2remove(i) = 1;
            elseif ~isempty(regexp(tline, 'Video onset', 'once'))
                [val, ~] = sscanf(tline, 'MSG %d Video onset %d');
                % remove catch trials for now, later maybe use (but I'll need to remove the test lines below I think)
                % alternative is to not use them and base the RDMs on only non-catch trials, like it is now I think
                if val(2) > 100
                    lines2remove(i) = 1;
                end 
            end           
        end
        
        aline(logical(lines2remove)) = [];
        
        allTrials = struct('onset',[],'cond',[],'idx',[]);
        itrial=0;
        for i=1:numel(aline)
            tline = aline{i};
                
            if regexp(tline, 'Video onset')
                itrial=itrial+1;
                
                [val, ~] = sscanf(tline, 'MSG %d Video onset %d');
                allTrials(itrial).onset = val(1);
                allTrials(itrial).cond = val(2);
                allTrials(itrial).idx = i;
            else
                % all other lines are not parsed
            end
            
        end
        
        PTBtriggers = extractfield(output,'trigger');
        PTBtriggers(PTBtriggers>100) = [];
        if max(PTBtriggers-extractfield(allTrials,'cond')) > 0% compare trial codes to PTB data
            error('PTB and eyetracker trial codes do not match!')
        end
        
        % fill a fieldtrip dataset
        ft_ds = [];
        checkIDtrials = zeros(length(allTrials),1);
        for itrial = 1:length(allTrials)
            ft_ds.sampleinfo(itrial,:) = [allTrials(itrial).onset-.2*fsample allTrials(itrial).onset+5*fsample];
            ft_ds.trialinfo(itrial,1) = allTrials(itrial).cond;
            samplecount=0;
            count = 0;
            for isample = allTrials(itrial).idx-.2*fsample:allTrials(itrial).idx+5*fsample

                tline = aline{isample};

                % remove video onset message
                if strcmp(tline(1:3),'MSG')
                    count = count+1;
                    continue
                end
                    
                samplecount=samplecount + 1;

                % small check to see if an eye was lost, and if so, which eye:
                check = split(tline);
                checkID = ones(7,1);
                for i = 1:7
                    if str2double(check{i}) == 0 || isnan(str2double(check{i}))
                        checkID(i) = 0;
                    end
                end
                checkID = logical(checkID);
                
                % mark in which trials an eye was lost:
                if sum(checkID) < 7
                    checkIDtrials(itrial) = 1;
                end
                
                % now extract eyetracker data and convert to double:
                temp = str2double(check);
                temp(isnan(temp)) = 0;

                ft_ds.trial(itrial,:,samplecount) = temp(1:7);
                
            end% sample loop
            
            % check whether last minus first sample is indeed 5 sec long
            if ft_ds.trial(itrial,1,end)+1 - ft_ds.trial(itrial,1,1) ~= 5200
               error(['Trial ' num2str(itrial) ' is not exactly 5.2 seconds long!']);
            end            
                        
        end% trial loop
        
        % remove sample numbers
        ft_ds.trial(:,1,:) = [];
        
        % How many samples was the eye lost per trial, and change values
        % there to NaN, so when downsampling-by-averaging below these are
        % ignored. After that they're interpolated for MVPA.
        missedSampleNum{iSub,iRun} = sum(sum(ft_ds.trial(logical(checkIDtrials),:,:)==0,2),3);
        ft_ds.trial(ft_ds.trial==0) = NaN;
        
        %% "round" time dimension (sometimes there seems to be a mismatch between runs)
        
        %% downsample
        A =  reshape(ft_ds.trial, size(ft_ds.trial,1), size(ft_ds.trial,2), size(ft_ds.trial,3)/(size(ft_ds.trial,3)*cfg.downsample/1000),[]);
        ft_ds.trial = squeeze(nanmean(A,3));
        ft_ds.trial(isnan(ft_ds.trial)) = 0;%change back to zero after using nanmean, because I'm not sure what classifier does with nan, better to simply use zero for those very few datapoints
        
        % remove trials that still have too many samples (10%) missing
        trials2remove = sum(sum(ft_ds.trial==0,2),3)./numel(ft_ds.trial(1,:,:)) > .1;
        missedTrialNum(iSub,iRun) = sum(trials2remove);
        ft_ds.trial(trials2remove,:,:) = [];
        ft_ds.trialinfo(trials2remove) = [];
        ft_ds.sampleinfo(trials2remove,:) = [];
        
        %% specify other ft fields
        ft_ds.time = -.2:1/cfg.downsample:5-1/cfg.downsample;
        ft_ds.label = {'left_x';'left_y';'left_pupil';'right_x';'right_y';'right_pupil'};
        ft_ds.dimord = 'rpt_chan_time';
        
        ds{iRun} = cosmo_meeg_dataset(ft_ds); % convert to cosmo for stacking the datasets
        
    end% run loop
    
    % join runs, re-add sensor information
    ds_all = cosmo_stack(ds);
    data = cosmo_map2meeg(ds_all); % ..and back to fieldtrip
    save(sprintf('%s%cSUB%02d',outdir, filesep, iSub), 'data','missedSampleNum','missedTrialNum', '-v7.3');
end


