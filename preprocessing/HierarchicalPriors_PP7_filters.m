% Filters
clearvars;

% Input files
rootdir = '\\XXX\HierarchicalPriors\data\MEG\brainstorm_database\Unpredict\data\';

% Input files
subjectKEY = {1 'MROS' ; 2 'MRSN' ; 3 'LANR' ; 4 'HLKR' ; 5 'LCBR' ; 6 'MRDG' ; 7 'SMHS' ; 8 'DBBU' ; 9	'PTCR' ; 10 'CRSA' ; 11 'CICE' ;...
    12 'CIBN' ; 13 'FARN' ; 14 'RGBR' ; 15 'LRSU' ; 16 'NNGG' ; 17 'POML' ; 18 'LRVL' ; 19 'SLBG' ; 20 'AAMN' ; 21 'RTKM' ; 22 'SNMN' ;...
    23 'GIMR' ; 24 'ATFC' ; 25 'MRVN' ; 26 'LCTY' ; 27 'JNFN' ; 28 'RBTM' ; 29 'SLML' ; 30 'ANDA' ; 31 'LLMS' ; 32 'AKAL' ; 33 'TTMD' ;...
    34 'ANDL' ; 35 'NCGR' ; 36 'RTDC' ; 37 'PTCS' ; 38 'GGSI' ; 39 'AGCI' ; 40 'ANVC' ; 41 'DNCC' ; 42 'CAOM' ; 43 'MRBU'};

SubjectNames = subjectKEY(:,2)';

for isub = 10
    
    sub = SubjectNames{isub};
    runfolders = dir([rootdir sub filesep '*@raw*']);
    sFiles = cell(1,length(runfolders));
    for irun=1:length(runfolders)
        runfile = dir([rootdir sub filesep runfolders(irun).name filesep '*raw*']);
        sFiles{1,irun} = [sub filesep runfolders(irun).name filesep runfile.name];
    end
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: Notch filter: 50Hz 100Hz
    sFiles = bst_process('CallProcess', 'process_notch', sFiles, [], ...
        'sensortypes', 'MEG, EEG', ...
        'freqlist',    [50, 100], ...
        'cutoffW',     0.5, ...
        'useold',      0, ...
        'read_all',    0);
    
    % Process: Resample: 500Hz
    sFiles = bst_process('CallProcess', 'process_resample', sFiles, [], ...
        'freq',     500, ...
        'read_all', 0);
    
    % Process: Power spectrum density (Welch)
    sFiles = bst_process('CallProcess', 'process_psd', sFiles, [], ...
        'timewindow',  [], ...
        'win_length',  1, ...
        'win_overlap', 50, ...
        'sensortypes', 'MEG', ...
        'win_std',     0, ...
        'edit',        struct(...
        'Comment',         'Power', ...
        'TimeBands',       [], ...
        'Freqs',           [], ...
        'ClusterFuncTime', 'none', ...
        'Measure',         'power', ...
        'Output',          'all', ...
        'SaveKernel',      0));
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);

    %% remove notch-filtered data to save space, we anyway keep the downsampled data after notch-filtering
    % find correct filenames
    runfolders = dir([rootdir sub filesep '*notch']);
    sFiles = cell(1,length(runfolders));
    for irun=1:length(runfolders)
        runfile = dir([rootdir sub filesep runfolders(irun).name filesep '*raw*']);
        sFiles{1,irun} = [sub filesep runfolders(irun).name filesep runfile.name];
    end
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: Delete folders
    sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
        'target', 2);  % Delete folders
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);
  
    
end

% Average PSD over runs for channel quality check
for isub = [3 5:7 10 12:24 27:31 33 34 36:40]
    
    sub = SubjectNames{isub};
    
    rundirs = dir(fullfile(rootdir,sub,'*resample'));
    
    % Input files
    sFiles = cell(1,6);
    for irun = 1:length(rundirs)
        file = dir(fullfile(rootdir,sub,rundirs(irun).name,'*timefreq*'));
        sFiles{1,irun} = fullfile(sub,rundirs(irun).name,file.name);
    end

    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: Average: Everything
    sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
        'avgtype',       1, ...  % Everything
        'avg_func',      1, ...  % Arithmetic average:  mean(x)
        'weighted',      0, ...
        'matchrows',     1, ...
        'iszerobad',     1);
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);
    
end
