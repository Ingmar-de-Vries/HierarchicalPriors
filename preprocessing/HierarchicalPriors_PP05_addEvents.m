% add events
clearvars

% Input files
rootdir = '\\XXX\HierarchicalPriors\data\MEG\brainstorm_database\Unpredict\data\';

% Input files
subjectKEY = {1 'MROS' ; 2 'MRSN' ; 3 'LANR' ; 4 'HLKR' ; 5 'LCBR' ; 6 'MRDG' ; 7 'SMHS' ; 8 'DBBU' ; 9	'PTCR' ; 10 'CRSA' ; 11 'CICE' ;...
    12 'CIBN' ; 13 'FARN' ; 14 'RGBR' ; 15 'LRSU' ; 16 'NNGG' ; 17 'POML' ; 18 'LRVL' ; 19 'SLBG' ; 20 'AAMN' ; 21 'RTKM' ; 22 'SNMN' ;...
    23 'GIMR' ; 24 'ATFC' ; 25 'MRVN' ; 26 'LCTY' ; 27 'JNFN' ; 28 'RBTM' ; 29 'SLML' ; 30 'ANDA' ; 31 'LLMS' ; 32 'AKAL' ; 33 'TTMD' ;...
    34 'ANDL' ; 35 'NCGR' ; 36 'RTDC' ; 37 'PTCS' ; 38 'GGSI' ; 39 'AGCI' ; 40 'ANVC' ; 41 'DNCC' ; 42 'CAOM' ; 43 'MRBU'};

SubjectNames = subjectKEY(:,2)';

for isub = [18:24 27:31 33 34 36:40]
    
    runfolders = dir([rootdir SubjectNames{isub} filesep '*@raw*']);
    sFiles = cell(1,length(runfolders));
    for irun=1:length(runfolders)
        runfile = dir([rootdir SubjectNames{isub} filesep runfolders(irun).name filesep '*raw*']);
        sFiles{1,irun} = [SubjectNames{isub} filesep runfolders(irun).name filesep runfile.name];
    end
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: Read from channel
    sFiles = bst_process('CallProcess', 'process_evt_read', sFiles, [], ...
        'stimchan',  'STI102', ...
        'trackmode', 1, ...  % Value: detect the changes of channel value
        'zero',      0);
    
    % Process: Merge events
    sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
        'evtnames', '1,2,3,4,5,6,7,8,9,10,11,12,13,14,21,22,23,24,25,26,27,28,29,30,31,32,33,34,41,42,43,44,45,46,47,48,49,50,51,52,53,54,101,102,103,104,105,106,107,108,109,110,111,112,113,114,121,122,123,124,125,126,127,128,129,130,131,132,133,134,141,142,143,144,145,146,147,148,149,150,151,152,153,154', ...
        'newname',  'video_onset');
    
    % Process: Merge events
    sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
        'evtnames', '201,202,203,204,205,206,207,208,209,210,211,212,213,214,221,222,223,224,225,226,227,228,229,230,231,232,233,234,241,242,243,244,245,246,247,248,249,250,251,252,253,254', ...
        'newname',  'test_onset');
    
    % Process: Merge events
    sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
        'evtnames', '61,68', ...
        'newname',  'response');
    
    % Process: Read from channel
    sFiles = bst_process('CallProcess', 'process_evt_read', sFiles, [], ...
        'stimchan',  'STI102', ...
        'trackmode', 1, ...  % Value: detect the changes of channel value
        'zero',      0);
    
    % Process: Delete events
    sFiles = bst_process('CallProcess', 'process_evt_delete', sFiles, [], ...
        'eventname', '61,68');
    
    EOGchan = 'EEG063';
    
    % Process: Detect blink with custom settings
    sFiles = bst_process('CallProcess', 'process_evt_detect', sFiles, [], ...
        'eventname',    'blink', ...
        'channelname',  EOGchan, ...
        'timewindow',   [], ...
        'bandpass',     [0.1, 15], ...
        'threshold',    2, ...
        'blanking',     0.5, ...
        'isnoisecheck', 1, ...
        'isclassify',   0);
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);
    
end% subject loop
