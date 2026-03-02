%% ===== CONVERT EPOCHED DATA FROM BRAINSTORM TO FIELDTRIP =====
% this script exports epoched data in fieldtrip compatible format

clearvars;

% set paths
indirPTB = '\\XXX\HierarchicalPriors\data\MEG\PTBoutput';
outdir = '\XXX\HierarchicalPriors\data\MEG\PreprocessedSensor';
addpath('\\XXX\HierarchicalPriors\toolboxes\fieldtrip-20191113');
addpath('\\XXX\HierarchicalPriors\code\preprocessing');
ft_defaults

subjectKEY = {1 'MROS' ; 2 'MRSN' ; 3 'LANR' ; 4 'HLKR' ; 5 'LCBR' ; 6 'MRDG' ; 7 'SMHS' ; 8 'DBBU' ; 9	'PTCR' ; 10 'CRSA' ; 11 'CICE' ;...
    12 'CIBN' ; 13 'FARN' ; 14 'RGBR' ; 15 'LRSU' ; 16 'NNGG' ; 17 'POML' ; 18 'LRVL' ; 19 'SLBG' ; 20 'AAMN' ; 21 'RTKM' ; 22 'SNMN' ;...
    23 'GIMR' ; 24 'ATFC' ; 25 'MRVN' ; 26 'LCTY' ; 27 'JNFN' ; 28 'RBTM' ; 29 'SLML' ; 30 'ANDA' ; 31 'LLMS' ; 32 'AKAL' ; 33 'TTMD' ;...
    34 'ANDL' ; 35 'NCGR' ; 36 'RTDC' ; 37 'PTCS' ; 38 'GGSI' ; 39 'AGCI' ; 40 'ANVC' ; 41 'DNCC' ; 42 'CAOM' ; 43 'MRBU'};

%PTB subject files
subfilz = dir(fullfile(indirPTB,'*block1*'));

% set some parameters
cfg = [];
cfg.realign2photodiode = 1;
cfg.FTformat = 2;%1 = Fieldtrip data format after ft_preprocessing, 2 = Fieldtrip data format after ft_timelockanalysis
%% ===== CREATE FIELDTRIP STRUCTURE =====
% get the subjs from the actual database
s = bst_get('ProtocolSubjects');

%% ===== SUBJECTS LOOP =====

for iSubj = 11
    
    % subject number
    subnum = find(strcmp(s.Subject(iSubj).Name,subjectKEY(:,2)));
    
    % Extract Structure of Subj to process
    sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',  s.Subject(iSubj).Name, ...
        'tag',          'video_onset');
    
    % If you have to go back to change a preprocessing step, you could manually remove noisy trials here so you don't have to rerun the detection:
%     badfile = dir(fullfile(outdir,['badtrials_' s.Subject(iSubj).Name '*']));
%     load(fullfile(outdir,badfile.name));

%     sFiles(allBadTrials) = [];   
    
    FTdata	= struct;
    
    %% === LOOP THROUGH TRIALS ===
    
    for iTrial = 1:length(sFiles)
        
        % Export the brainstorm data structure of each trial into fieldtrip timelocked format
        [data, dataBST, ChannelMat] = out_fieldtrip_data(sFiles(iTrial).FileName, sFiles(iTrial).ChannelFile, ...
            'MEG,MISC008,STI102,EEG BAD LOC', 0);
        
        % Combining trials
        FTdata.time{1,iTrial}               = cell2mat(data.time);
        FTdata.trial{1,iTrial}				= cell2mat(data.trial);
        
        % Find trigger value
        for ievent = 1:length(dataBST.Events)
            if dataBST.Events(ievent).times == 0
                FTdata.trialinfo(iTrial,1) = str2double(dataBST.Events(ievent).label);
            end
        end

        % Some extra info not necessary for FT but handy nonetheless
        FTdata.BSTevents{iTrial}			= dataBST.Events;
        FTdata.original_trialnum{iTrial}             = dataBST.Comment;
        
    end
    
    % complete the FieldTrip Structure with common fields across trials
    FTdata.label							= data.label;
    FTdata.grad                             = data.grad;
    
    %% realign data to photodiode
    if cfg.realign2photodiode == 1 && ~any(subnum==[17,18])% photodiode not recorded for those two subjects
        FTdata = Unpredict_PP10_realign2photodiode(FTdata);
    end
    
    %% transform from FT preprocessed format to FT timelocked format, and remove some unnecessary fields
    % store the extra BST info because it gets lost when transforming to FT timelocked format
    BSTevents = FTdata.BSTevents;
    
    if cfg.FTformat == 2
        
        FTdata = rmfield(FTdata,{'BSTevents','original_trialnum'});

        ft_cfg=[];
        ft_cfg.keeptrials = 'yes';
        ft_cfg.removemean = 'no';
        ft_cfg.channel = 'meg';
        FTdata = ft_timelockanalysis(ft_cfg,FTdata);
    end
    
    FTdata.BSTevents = BSTevents;
    
    %% check triggers with PTB file and fix if a trigger value is strange  
    runfilz = dir(sprintf('%s%cUnpredict_subject%02d_block*',indirPTB, filesep, subnum));
    
    data = [];
    for irun = 1:length(runfilz)
    
        load(fullfile(indirPTB,runfilz(irun).name));

        if irun == 1
        	data = output;
        else
            data = [data output];
        end
        
    end
%     data(allBadTrials) = [];
    
    PTBtriggers=extractfield(data,'trigger')';
    MEGtriggers = FTdata.trialinfo;
    
    trigger_mismatch = find(PTBtriggers~=MEGtriggers);        
%     FTdata.trialinfo(trigger_mismatch) = PTBtriggers(trigger_mismatch);
    
    %% ===== WRITE THE SUBJECT STRUCTURE TO DISK =====
    save(sprintf('%s%cpreprocessedSensor_SUB%02d%s.mat', outdir, filesep, subnum, s.Subject(iSubj).Name), 'FTdata', '-v7.3');
    save(['XXX\HierarchicalPriors\progress_sub' num2str(iSubj)],'cfg');
    
    % if mismatch between PTB and MEG triggers, manually check for this subject why that might be the case, but continue loop with others first
    if ~isempty(trigger_mismatch)
        save([outdir filesep 'trigger_mismatch_' s.Subject(iSubj).Name, '.mat'], 'trigger_mismatch','PTBtriggers','MEGtriggers','-v7.3');
    end
    
end
