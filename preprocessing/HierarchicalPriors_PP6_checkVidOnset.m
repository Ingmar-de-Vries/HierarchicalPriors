
% In some runs we find more than 84 video onset events (or 108 if including catch trials), because of duplicates.
% In this piece of code we search for those duplicates. 
% To double check you can compare it to the trigger values in the PTB output.
% After identifying them, remove them manually in the BST gui.
% Ingmar de Vries October 2020
clearvars;

% set paths
addpath('\\XXX\HierarchicalPriors\code\preprocessing');
filedir = '\\XXX\HierarchicalPriors\data\MEG\manualEventFixes\';
file = dir(fullfile(filedir, '*mrk*'));

file = file(5)

% Load BST events
BSTevents = readBSTevents([filedir filesep file.name]);
BSTevents(isnan(BSTevents)) = [];
BSTevents = sort(BSTevents);

for itrial = 2:length(BSTevents)
    if BSTevents(itrial) <= BSTevents(itrial-1)+5
        disp(['Dulicate is trial ' num2str(itrial) ' at t = ' num2str(BSTevents(itrial))]);
    end
end