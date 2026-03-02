%% ===== APPLY MINIMUM NORM ESTIMATION AND SAVE RESULTING INVERSION KERNEL =====
% This script applies minimum norm estimation (compute head model + noise
% covariance + data convariance --> MNE), and stores the resulting
% inversion kernel

% AUTHOR: Ingmar de Vries, September 2021
clearvars;

% output directory
outdir = '\\XXX\HierarchicalPriors\data\MEG\InversionKernels';
addpath('\\XXX\HierarchicalPriors\code\preprocessing');
badtrialdir = '\\XXX\HierarchicalPriors\data\MEG\PreprocessedSensor';
indirPTB = '\\XXX\HierarchicalPriors\data\MEG\PTBoutput';
indirBST = '\\XXX\HierarchicalPriors\data\MEG\brainstorm_database\Unpredict\data';

subjectKEY = {1 'MROS' ; 2 'MRSN' ; 3 'LANR' ; 4 'HLKR' ; 5 'LCBR' ; 6 'MRDG' ; 7 'SMHS' ; 8 'DBBU' ; 9	'PTCR' ; 10 'CRSA' ; 11 'CICE' ;...
    12 'CIBN' ; 13 'FARN' ; 14 'RGBR' ; 15 'LRSU' ; 16 'NNGG' ; 17 'POML' ; 18 'LRVL' ; 19 'SLBG' ; 20 'AAMN' ; 21 'RTKM' ; 22 'SNMN' ;...
    23 'GIMR' ; 24 'ATFC' ; 25 'MRVN' ; 26 'LCTY' ; 27 'JNFN' ; 28 'RBTM' ; 29 'SLML' ; 30 'ANDA' ; 31 'LLMS' ; 32 'AKAL' ; 33 'TTMD' ;...
    34 'ANDL' ; 35 'NCGR' ; 36 'RTDC' ; 37 'PTCS' ; 38 'GGSI' ; 39 'AGCI' ; 40 'ANVC' ; 41 'DNCC' ; 42 'CAOM' ; 43 'MRBU'};

% get the subjs from the brainstorm database
s = bst_get('ProtocolSubjects');

% remove group level data
s.Subject(contains(extractfield(s.Subject,'Name'),'Group')) = [];

%% ===== SUBJECTS LOOP =====
trialnum = zeros(length(s.Subject),1);
for iSubj = 1:numel(s.Subject)
    
    % Extract Structure of Subj to process
    sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',  s.Subject(iSubj).Name, ...
        'tag',          'video_onset');
        
    % select trials not to include in source reconstruction (i.e., trials with noise during the baseline interval which is used for the noice
    % covariance matrix)
    clear allBadTrials
    load(fullfile(badtrialdir,['badtrials_' s.Subject(iSubj).Name '.mat']),'allBadTrials');
    allBadTrials = find(allBadTrials);
    
    trials2remove = false(length(sFiles),1);
    for ibadtrial = 1:length(allBadTrials)
        
        % load bad trial info
        load(fullfile(indirBST,sFiles(allBadTrials(ibadtrial)).FileName),'Events');
        
        badID = contains(extractfield(Events,'label'),'bad');
        badID = find(badID);
        
        for ibadtype = 1:length(badID)
            
            bad_times = Events(badID(ibadtype)).times;
            
            for ibadseg = 1:size(bad_times,2)
                
                % we want to remove noisy trials in 3 cases:
                % if start of bad segment falls within the baseline window (-0.5 to 0)
                % if end of bad segment falls within the baseline interval
                % OR if start of bad segment falls before the interval, AND end falls after the interval, i.e., the bad segment encompasses the
                % whole baseline interval
                if bad_times(1,ibadseg) > -0.5 && bad_times(1,ibadseg) < 0 || bad_times(2,ibadseg) > -0.5 && bad_times(2,ibadseg) < 0 || bad_times(1,ibadseg) < -0.5 && bad_times(2,ibadseg) > 0
                    trials2remove(allBadTrials(ibadtrial)) = true;
                end
                
            end
            
        end
        
    end
    
    % now remove files
    sFiles(trials2remove) = [];
    trialnum(iSubj) = length(sFiles);
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: Compute head model
    sFiles = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
        'Comment',     '', ...
        'sourcespace', 1, ...  % Cortex surface
        'meg',         3, ...  % Overlapping spheres
        'eeg',         3, ...  % OpenMEEG BEM
        'ecog',        2, ...  % OpenMEEG BEM
        'seeg',        2, ...  % OpenMEEG BEM
        'openmeeg',    struct(...
        'BemFiles',     {{}}, ...
        'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
        'BemCond',      [1, 0.0125, 1], ...
        'BemSelect',    [1, 1, 1], ...
        'isAdjoint',    0, ...
        'isAdaptative', 1, ...
        'isSplit',      0, ...
        'SplitLength',  4000), ...
        'channelfile', '');
    
    % Process: Compute covariance (noise or data)
    sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
        'baseline',       [-1, -0.002], ...
        'datatimewindow', [0, 0], ...
        'sensortypes',    'MEG', ... %'MEG GRAD', ...
        'target',         1, ...  % Noise covariance     (covariance over baseline time window)
        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'copymatch',      0, ...
        'replacefile',    1);  % Replace
    
    % Process: Compute sources [2018]
    sFiles = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
        'output',  1, ...  % Kernel only: shared
        'inverse', struct(...
        'Comment',        'MN: MEG ALL', ... %'sLORETA: MEG ALL', ... %'MN: MEG GRAD', ...
        'InverseMethod',  'minnorm', ...
        'InverseMeasure', 'amplitude', ... %'sloreta', ... %
        'SourceOrient',   {{'fixed'}}, ... % {{'free'}}, ...
        'Loose',          0.2, ...
        'UseDepth',       0, ... %1, ...
        'WeightExp',      0.5, ...
        'WeightLimit',    10, ...
        'NoiseMethod',    'reg', ... %'median', ... %
        'NoiseReg',       0.1, ...
        'SnrMethod',      'fixed', ...
        'SnrRms',         1e-06, ...
        'SnrFixed',       3, ...
        'ComputeKernel',  1, ...
        'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));%{{'MEG GRAD'}}));
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);
    
end