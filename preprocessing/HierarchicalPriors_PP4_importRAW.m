% import (link to) raw MEG data
clearvars

rootdir = '\\XXX\HierarchicalPriors\data\MEG\ELEKTAoutput';

% Input files
subjectKEY = {1 'MROS' ; 2 'MRSN' ; 3 'LANR' ; 4 'HLKR' ; 5 'LCBR' ; 6 'MRDG' ; 7 'SMHS' ; 8 'DBBU' ; 9	'PTCR' ; 10 'CRSA' ; 11 'CICE' ;...
    12 'CIBN' ; 13 'FARN' ; 14 'RGBR' ; 15 'LRSU' ; 16 'NNGG' ; 17 'POML' ; 18 'LRVL' ; 19 'SLBG' ; 20 'AAMN' ; 21 'RTKM' ; 22 'SNMN' ;...
    23 'GIMR' ; 24 'ATFC' ; 25 'MRVN' ; 26 'LCTY' ; 27 'JNFN' ; 28 'RBTM' ; 29 'SLML' ; 30 'ANDA' ; 31 'LLMS' ; 32 'AKAL' ; 33 'TTMD' ;...
    34 'ANDL' ; 35 'NCGR' ; 36 'RTDC' ; 37 'PTCS' ; 38 'GGSI' ; 39 'AGCI' ; 40 'ANVC' ; 41 'DNCC' ; 42 'CAOM' ; 43 'MRBU'};

SubjectNames = subjectKEY(:,2)';

for isub = [3 5 7 10 12 13 14 15 17 20 21 28 31 33 34 36 37 39 40]
    
    rundir = fullfile(rootdir,SubjectNames{isub});
    filz = dir(fullfile(rundir,'*trans_tsss.fif'));
    
    sFiles = [];
    
    RawFiles = cell(1,length(filz));
    for ifile = 1:length(filz)
        RawFiles{ifile} = fullfile(rundir,filz(ifile).name);
    end
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: Create link to raw file
    sFiles = bst_process('CallProcess', 'process_import_data_raw', sFiles, [], ...
        'subjectname',    SubjectNames{isub}, ...
        'datafile',       {RawFiles, 'FIF'}, ...
        'channelreplace', 0, ...
        'channelalign',   1, ...
        'evtmode',        'value');
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);

end
