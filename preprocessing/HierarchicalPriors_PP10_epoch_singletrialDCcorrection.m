%% import data into database and epoch
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
    
    runfolders = dir([rootdir sub filesep '*resample']);
    runfolders = runfolders(contains(extractfield(runfolders,'name'),'@raw'));
    sFiles = cell(1,length(runfolders));
    for irun=1:length(runfolders)
        runfile = dir([rootdir sub filesep runfolders(irun).name filesep '*raw*']);
        sFiles{1,irun} = [sub filesep runfolders(irun).name filesep runfile.name];
    end
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: Import MEG/EEG: Events
    sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
        'subjectname', sub, ...
        'condition',   '', ...
        'eventname',   'video_onset', ...%
        'timewindow',  [], ...
        'epochtime',   [-1.5, 6.5], ...
        'createcond',  0, ...
        'ignoreshort', 1, ...
        'usectfcomp',  0, ...
        'usessp',      1, ...
        'freq',        [], ...
        'baseline',    []);
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);
    
end

%% single-trial baseline correction (DC offset)
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
    
    runfolders = dir([rootdir sub filesep '*resample']);
    folders2keep = [];
    for irun = 1:length(runfolders)
        if ~contains(runfolders(irun).name,'@raw')
            folders2keep = [folders2keep irun];
        end
    end
    runfolders = runfolders(folders2keep);
    
    sFiles = cell(1,108*length(runfolders));
    for irun=1:length(runfolders)
        filz = dir([rootdir sub filesep runfolders(irun).name filesep '*trial*']);
        
        for ifile = 1:length(filz)
            sFiles{1,(irun-1)*length(filz)+ifile} = [sub filesep runfolders(irun).name filesep filz(ifile).name];
        end
    end
    
    sFiles = sFiles(~cellfun('isempty', sFiles));%if less than the initialized 648 trials, remove empty cells here
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: DC offset correction: [-500ms,-2ms]
    sFiles = bst_process('CallProcess', 'process_baseline', sFiles, [], ...
        'baseline',    [-0.5, -0.002], ...
        'sensortypes', 'MEG', ...
        'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
        'overwrite',   1);
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);
 
end

%% bad trials
% visually inspect whether the automatically detected trials should indeed all be removed, 
% e.g., sometimes bad segment is in our padding windows, or might not be bad
% in our eyes, so mark those below
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
    
    runfolders = dir([rootdir sub filesep '*resample']);
    folders2keep = [];
    for irun = 1:length(runfolders)
        if ~contains(runfolders(irun).name,'@raw')
            folders2keep = [folders2keep irun];
        end
    end
    runfolders = runfolders(folders2keep);

    trials2reject = false(108,length(runfolders));
    for irun = 1:length(runfolders)
       
        trials = dir([rootdir sub filesep runfolders(irun).name filesep 'data_video_onset_trial*']);
        for itrial = 1:length(trials)
            
            % load event info for single trial
            load([rootdir sub filesep runfolders(irun).name filesep trials(itrial).name],'Events');
            
            % for catch trials, set end of trial at test onset
            % for normal trials set end of trial at video offset, i.e., 5 sec
            if any(strcmp(extractfield(Events,'label'),'test_onset'))
                trialend = extractfield(Events(strcmp(extractfield(Events,'label'),'test_onset')),'times');
                trialend = trialend(1);
            else
                trialend = 5;
            end
            
            % BST applies some padding, but not clear from documentation how much. 
            % We actually don't care that the padding falls in our long [-1 5] interval, 
            % but only the actually noisy part. So we remove the padding again: 
            padding = 0.3;
            trialstart = -1+padding;
            trialend = trialend-padding;
            
            % check if current trial has a bad segment
            badseg = find(contains(extractfield(Events,'label'),'bad')); 
            
            % loop over all types bad segments, i.e., 1-7Hz or 40-240Hz are stored separately
            for iseg = 1:length(badseg)
                
                % find timing of bad segments
                badtimes = Events(badseg(iseg)).times;
                for ibad = 1:size(badtimes,2)
                    % we want to keep the bad segment marking in OR 3 cases:
                    % if start of bad segment falls within our interval of interest
                    % if end of bad segment falls within our interval of interest
                    % OR if start of bad segment falls before the interval, AND end falls after the interval, i.e., the bad segment encompasses the
                    % whole trial
                    if (badtimes(1,ibad) > trialstart && badtimes(1,ibad) < trialend) || (badtimes(2,ibad) > trialstart && badtimes(2,ibad) < trialend) || (badtimes(1,ibad) < trialstart && badtimes(2,ibad) > trialend)
                        trials2reject(itrial,irun) = true;
                    end
                end
            
            end
        end

    end
end

% visual inspection in the BST GUI, and comparison with the remaining list of bad trials displayed in command window:
irun = 6;
find(trials2reject(:,irun))'

% if we want to deselect additional bad segments after visual inspection, manually list them here
seg2keep = [15 16 18 28 44 45 47 50 62 69 70 71 92 94 104];
trials2reject(seg2keep,irun) = 0;

% now save the trials2reject after visual inspection
allBadTrials = reshape(trials2reject,[],1);
badtrialdir = '\\XXX\HierarchicalPriors\data\MEG\PreprocessedSensor';
save(fullfile(badtrialdir,['badtrials_' sub]),'allBadTrials');

%% In GUI: select all imported data folders of all runs, right click, group folders.

% %% subject CRSA 
% seg2keep = [43 60 80 87 90 92 100 105:108];
% seg2keep = [2 4 17 24 34 37 39 41 48 63 65 70 82 101 106 108]; 
% seg2keep = [13 46 63 67 68 72 76 77 89 90];
% seg2keep = [1 6 16 18 22 23 26 27 33 46 47 53 54 68 70 72 79 82 86 87 98 100 103];
% seg2keep = [15 16 22 34 57 71 72 84 89 93 98 100 101 105 106 108];
% seg2keep = [15 16 18 28 44 45 47 50 62 69 70 71 92 94 104];

