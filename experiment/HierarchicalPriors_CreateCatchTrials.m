%% Create pool of all potential catch trials
CatchTrials = cell(4,28);
% first dimension (4 rows) = 
%   1. sequence number (14x)
%   2. stop time (sec)
%   3. sequence after occlusion
%   4. correct answer (1 = correct, 0 = incorrect) 

%% kinematics
% sequence 1
CatchTrials{1,1} = 1;
CatchTrials{2,1} = 3.3;
CatchTrials{3,1} = 1;
CatchTrials{4,1} = 1;

CatchTrials{1,2} = 1;
CatchTrials{2,2} = 3.3;
CatchTrials{3,2} = 3;
CatchTrials{4,2} = 0;

% sequence 2
CatchTrials{1,3} = 2;
CatchTrials{2,3} = 4.2;
CatchTrials{3,3} = 2;
CatchTrials{4,3} = 1;

CatchTrials{1,4} = 2;
CatchTrials{2,4} = 4.2;
CatchTrials{3,4} = 3;
CatchTrials{4,4} = 0;

% sequence 3
CatchTrials{1,5} = 3;
CatchTrials{2,5} = 3.5;
CatchTrials{3,5} = 3;
CatchTrials{4,5} = 1;

CatchTrials{1,6} = 3;
CatchTrials{2,6} = 3.5;
CatchTrials{3,6} = 1;
CatchTrials{4,6} = 0;

% sequence 4
CatchTrials{1,7} = 4;
CatchTrials{2,7} = 4.4;
CatchTrials{3,7} = 4;
CatchTrials{4,7} = 1;

CatchTrials{1,8} = 4;
CatchTrials{2,8} = 4.4;
CatchTrials{3,8} = 3;
CatchTrials{4,8} = 0;

% sequence 5
CatchTrials{1,9} = 5;
CatchTrials{2,9} = 3.2;
CatchTrials{3,9} = 5;
CatchTrials{4,9} = 1;

CatchTrials{1,10} = 5;
CatchTrials{2,10} = 3.2;
CatchTrials{3,10} = 6;
CatchTrials{4,10} = 0;

% sequence 6
CatchTrials{1,11} = 6;
CatchTrials{2,11} = 4.2;
CatchTrials{3,11} = 6;
CatchTrials{4,11} = 1;

CatchTrials{1,12} = 6;
CatchTrials{2,12} = 4.2;
CatchTrials{3,12} = 5;
CatchTrials{4,12} = 0;

% sequence 7
CatchTrials{1,13} = 7;
CatchTrials{2,13} = 3.4;
CatchTrials{3,13} = 7;
CatchTrials{4,13} = 1;

CatchTrials{1,14} = 7;
CatchTrials{2,14} = 3.4;
CatchTrials{3,14} = 8;
CatchTrials{4,14} = 0;

% sequence 8
CatchTrials{1,15} = 8;
CatchTrials{2,15} = 3.4;
CatchTrials{3,15} = 8;
CatchTrials{4,15} = 1;

CatchTrials{1,16} = 8;
CatchTrials{2,16} = 3.4;
CatchTrials{3,16} = 7;
CatchTrials{4,16} = 0;

% sequence 9
CatchTrials{1,17} = 9;
CatchTrials{2,17} = 4.2;
CatchTrials{3,17} = 9;
CatchTrials{4,17} = 1;

CatchTrials{1,18} = 9;
CatchTrials{2,18} = 4.2;
CatchTrials{3,18} = 7;
CatchTrials{4,18} = 0;

% sequence 10
CatchTrials{1,19} = 10;
CatchTrials{2,19} = 3.3;
CatchTrials{3,19} = 10;
CatchTrials{4,19} = 1;

CatchTrials{1,20} = 10;
CatchTrials{2,20} = 3.3;
CatchTrials{3,20} = 11;
CatchTrials{4,20} = 0;

% sequence 11
CatchTrials{1,21} = 11;
CatchTrials{2,21} = 3.3;
CatchTrials{3,21} = 11;
CatchTrials{4,21} = 1;

CatchTrials{1,22} = 11;
CatchTrials{2,22} = 3.3;
CatchTrials{3,22} = 10;
CatchTrials{4,22} = 0;

% sequence 12
CatchTrials{1,23} = 12;
CatchTrials{2,23} = 3.4;
CatchTrials{3,23} = 12;
CatchTrials{4,23} = 1;

CatchTrials{1,24} = 12;
CatchTrials{2,24} = 3.4;
CatchTrials{3,24} = 11;
CatchTrials{4,24} = 0;

% sequence 13
CatchTrials{1,25} = 13;
CatchTrials{2,25} = 2.9;
CatchTrials{3,25} = 13;
CatchTrials{4,25} = 1;

CatchTrials{1,26} = 13;
CatchTrials{2,26} = 2.9;
CatchTrials{3,26} = 12;
CatchTrials{4,26} = 0;

% sequence 14
CatchTrials{1,27} = 14;
CatchTrials{2,27} = 3.1;
CatchTrials{3,27} = 14;
CatchTrials{4,27} = 1;

CatchTrials{1,28} = 14;
CatchTrials{2,28} = 3.1;
CatchTrials{3,28} = 13;
CatchTrials{4,28} = 0;


save('C:\temporary_experiment\experiment\HierarchicalPriors_CatchTrialPool','CatchTrials');




