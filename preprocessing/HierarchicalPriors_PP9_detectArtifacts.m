% Detect additional artifacts: 1-7 Hz high amplitude for movements
% (incl. eye movements), 40-240 Hz high amplitude for muscle activity
clearvars; 

% Input files
rootdir = '\\\XXX\HierarchicalPriors\data\MEG\brainstorm_database\Unpredict\data\';

% Input files
subjectKEY = {1 'MROS' ; 2 'MRSN' ; 3 'LANR' ; 4 'HLKR' ; 5 'LCBR' ; 6 'MRDG' ; 7 'SMHS' ; 8 'DBBU' ; 9	'PTCR' ; 10 'CRSA' ; 11 'CICE' ;...
    12 'CIBN' ; 13 'FARN' ; 14 'RGBR' ; 15 'LRSU' ; 16 'NNGG' ; 17 'POML' ; 18 'LRVL' ; 19 'SLBG' ; 20 'AAMN' ; 21 'RTKM' ; 22 'SNMN' ;...
    23 'GIMR' ; 24 'ATFC' ; 25 'MRVN' ; 26 'LCTY' ; 27 'JNFN' ; 28 'RBTM' ; 29 'SLML' ; 30 'ANDA' ; 31 'LLMS' ; 32 'AKAL' ; 33 'TTMD' ;...
    34 'ANDL' ; 35 'NCGR' ; 36 'RTDC' ; 37 'PTCS' ; 38 'GGSI' ; 39 'AGCI' ; 40 'ANVC' ; 41 'DNCC' ; 42 'CAOM' ; 43 'MRBU'};

SubjectNames = subjectKEY(:,2)';

for isub = 10
    
    sub = SubjectNames{isub};
    
    runfolders = dir([rootdir sub filesep '*resample*']);
    sFiles = cell(1,length(runfolders));
    for irun=1:length(runfolders)
        runfile = dir([rootdir sub filesep runfolders(irun).name filesep '*raw*']);
        sFiles{1,irun} = [sub filesep runfolders(irun).name filesep runfile.name];
    end
    
    % Start a new report
    bst_report('Start', sFiles);

    % Process: Detect low freq noise, i.e., blinks, movements, etc.
    sFiles = bst_process('CallProcess', 'process_evt_detect_badsegment', sFiles, [], ...
        'timewindow',  [], ...
        'sensortypes', 'MEG', ...
        'threshold',   5, ... % (we don't actually care too much about these unless they're really large, so don't go too low with threshold)
        'isLowFreq',   1, ...
        'isHighFreq',  0);

        % Process: Detect high freq noise, i.e., muscle activity
    sFiles = bst_process('CallProcess', 'process_evt_detect_badsegment', sFiles, [], ...
        'timewindow',  [], ...
        'sensortypes', 'MEG', ...
        'threshold',   5, ... % (here we care a bit more about and generally muscle activity is estimated quite well 
        'isLowFreq',   0, ...
        'isHighFreq',  1);
    
%     Process: Rename events
    sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
        'src',   '40-240Hz', ...
        'dest',  'bad_40-240Hz');
    
    % Process: Rename events
    sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
        'src',   '1-7Hz', ...
        'dest',  'bad_1-7Hz');
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);

end
