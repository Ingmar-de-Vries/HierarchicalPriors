% Clear the workspace
close all;
clear all;
clc;
sca;

% For debugging purposes this makes experiment transparent. Simply comment
% for real experiment.
debug = false;
if debug
    
    PsychDebugWindowConfiguration([],.4);
    Screen('Preference', 'SkipSyncTests', 1);%don't do screen sync + size tests
end
Screen('Preference', 'SkipSyncTests', 0);%do screen sync + size tests
%----------------------------------------------------------------------
%                        INITIALIZE EXPERIMENT
%----------------------------------------------------------------------

% SOME GENERAL SETTINGS
rootdir = 'rootdirectory';
addpath(genpath(rootdir));
trackeye = false;
MEG = true;
if MEG
    trackeye = true;
end
if trackeye
    Screen('Preference', 'SkipSyncTests', 0);
end

% Number of repetitions of unique trials/sequences * conditions per block (14 unique sequences * 3 conditions [normal, upside, scrambled])
trialsPerCombination = 2;

% determine amount of practice trials.
trialAmountPractice = 60;

% FILE NAME, ETC.
rng('default');
rng('shuffle'); %seed random number generator

experiment_name = 'HierarchicalPriors';
subject_num = input('Please indicate subject number: ');
practice = input('Is this a practice run? [1/0]');
if practice
    experiment_name = [experiment_name '_practice'];
    block = 0;
else
    block = input('Please indicate number of real block: ');
end

if subject_num < 10
    filename = [experiment_name '_subject0' num2str(subject_num) '_block' num2str(block)];
else
    filename = [experiment_name '_subject' num2str(subject_num) '_block' num2str(block)];
end

savedir = [rootdir 'data\'];
addpath([rootdir 'experiment']);

if ~practice && exist([savedir filename '.mat'],'file')
    warning('Filename already exists, please make sure you have indicated the correct subject and block number!');
    return
end

% Monitor refresh rate, in MEG lab 120, in behavioural lab probably 60 but check
% Also, in MEG lab we may want to set it to a multiple of 50 Hz because of the videos (i.e. with 100 Hz, each video frame stays for 2 monitor frames)
% refrate = 59; % now we determine the refrate below with the PTB function: Screen('GetFlipInterval',window);
display.dist = 58;% distance from screen (cm)
display.width = 41; %width of screen (cm)
display.resolution = 1280; % number of pixels of display in horizontal direction

% PSYCHTOOLBOX
PsychDefaultSetup(2);
screens = Screen('Screens');
screenNumber = 1;%max(screens);

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseVirtualFramebuffer');

% define colors
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
bg = grey*0.8;
colors = [27 158 119;217 95 2;117 112 179;231,41,138;102,166,30]/255;%RGB values of colors that are colorblind and grayscale friendly
colorPlaceHolder = [0 0 0];

% Open the experiment window
[window, windowRect] = PsychImaging('OpenWindow',  screenNumber, bg, [], 32, 2);
refrate = 1/Screen('GetFlipInterval', window);
refrate = round(refrate);
% refrate = 40;

% Set PTB to top priority so no other running processes on this PC interfer
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

% Flip to clear
Screen('Flip', window);

% Get screen info
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[screenXpixels, screenYpixels] = Screen('WindowSize', window);%get screen size
[xCenter, yCenter] = RectCenter(windowRect);% get screen center

% For real experiment hide cursor
if ~debug
    HideCursor;
end

% STIMULI & CONDITION MATRIX
% define stimuli sizes, screen locations, and load videos
Screen('TextSize', window, 28);
fixCrossDimPix = angle2pix(display,.3);
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 2;

%sizes in visual angle divided by deg_per_pix to convert to pixels on screen
stimHeight= 400;%angle2pix(display,3.3);
stimWidth = 376;
rect4stim = [0 0 stimWidth stimHeight];
%stimPos 1 = left, 2 = right, 3 = center
stimPos = CenterRectOnPointd(rect4stim, xCenter, yCenter);

%load all videos here
stimdir = [rootdir 'experiment\Stimuli\videos'];
videos = dir(fullfile(stimdir, '*vid_120Hz.mp4'));
if subject_num < 10 
    fn_scramindices = fullfile(savedir,['HierarchicalPriors_subject0' num2str(subject_num) '_scramindices.mat']);
    fn_backup = fullfile(savedir,'scramidx_backup',['HierarchicalPriors_subject0' num2str(subject_num) '_scramindices_final.mat']);
else
    fn_scramindices = fullfile(savedir,['HierarchicalPriors_subject' num2str(subject_num) '_scramindices.mat']);
    fn_backup = fullfile(savedir,'scramidx_backup',['HierarchicalPriors_subject' num2str(subject_num) '_scramindices_final.mat']);
end
stimnum = length(videos);
framenum = 5*refrate;
stimuli = cell(3,stimnum,framenum);% cell with size: [conditions x sequences x frames]

DrawFormattedText(window,...
    'Please try to relax while the videos are loading...',...
    'center', yCenter/3, black);
Screen('Flip', window);

samples_all = cell(stimnum,1);
randID_all = cell(stimnum,1);
for istim = 1:stimnum
    
    %load videos and change frames to textures
    videoName = fullfile(stimdir,videos(istim).name);
    videoHeader = VideoReader(videoName);
    video = read(videoHeader);
    
    % if first run, create scramble indices for this subject and save
    % if second run or later, load the previously created indices
    if practice
        
        % piecewise scrambled, variable uniformly distributed segment lengths between:        
        minframes = round(0.2*videoHeader.framerate);% 200 msec = 24 frames | 240 msec = 
        maxframes = round(0.5*videoHeader.framerate);% 500 msec = 60 frames | 740 msec = 
        segsamp = round(minframes + (maxframes-minframes).*rand(100,1));
        segsamp = segsamp(1:find(cumsum(segsamp)>size(video,4),1)-1);
        segsamp(end+1) = size(video,4)-sum(segsamp);% last one will be to fill up rest of the video
        if sum(segsamp(end-1:end)) < maxframes% if last two segments together are still shorter than maxframes, just add them to prevent having a very short one
            segsamp(end-1) = sum(segsamp(end-1:end));
            segsamp(end) = [];
        end
        segsamp = [0 ; cumsum(segsamp)];
        samples = [segsamp(1:end-1)+1 segsamp(2:end)];
        
        randID = randperm(length(segsamp)-1);% random indexing, first try
        while any(diff(randID)==1)% if any subsequent segments are still in the original order, retry randomization, until it's not, so we don't accidentally have a much longer segment
            randID = randperm(length(segsamp)-1);
        end
        
        samples_all{istim} = samples;
        randID_all{istim} = randID;
        save(fn_scramindices,'samples_all','randID_all');
    else
        load(fn_scramindices);
        samples = samples_all{istim};
        randID = randID_all{istim};
    end

    % then do scrambling of video:
    newidx = [];
    for isegment = 1:length(samples)
        newidx = [newidx samples(randID(isegment),1):samples(randID(isegment),2)];
    end    
       
    for iframe = 1:framenum
        stimuli{1,istim,iframe} = Screen('MakeTexture', window, squeeze(video(:,:,:,iframe)));
        stimuli{2,istim,iframe} = Screen('MakeTexture', window, squeeze(video(end:-1:1,:,:,iframe)));
    end
    clear video
    stimuli(3,istim,:) = stimuli(1,istim,newidx);
end

% All variables
% video sequence
sequence = 1:stimnum;
conditions = 1:3;

% Make the matrix which will determine our condition combinations
% Row 1 = video sequence, row 2 = condition, row 3 = catch or no catch trial
condMatrixBase = repmat(sequence, 1, length(conditions));
condMatrixBase(2,:) = sort(repmat(conditions, 1, length(sequence)));

% Duplicate the condition matrix to get the full number of trials
condMatrix = repmat(condMatrixBase, 1, trialsPerCombination);
condMatrix(3,:) = 0;% no catch trials yet

% Add catch or no catch trial variable
% of the total amount of catch trials in this block: 16 = occlusion task, 8 = fixation cross task, i.e., 84 normal trials 
% catchtrial: 1 = catch trial, 2 = normal trial
if practice
    numcatch = [28 16];% occlusion trials and fixation trials
else
    numcatch = [16 8];%[16 8];
end

load([rootdir 'experiment\HierarchicalPriors_CatchTrialPool'],'CatchTrials');
% Times in CatchTrials variable indicate roughly middle of larger motion, so occlusion can start a little earlier and end a little later:
occlusionStart = 0.2;% seconds before manually picked middle of a larger motion
CatchTrials(2,:) = num2cell(cell2mat(CatchTrials(2,:)) - occlusionStart);

catchTrialPool = CatchTrials(:,randperm(size(CatchTrials,2),numcatch(1)));
for i=numcatch(1)+1:sum(numcatch)
    catchTrialPool{1,i} = randperm(stimnum,1);
    catchTrialPool{2,i} = 0.4 + 4.4 * rand(1,1);
    catchTrialPool{3,i} = 0;
    catchTrialPool{4,i} = 1;% correct response is always left for fixation cross task
end
numcatch = size(catchTrialPool,2);

catchTrialPool = catchTrialPool(:,randperm(numcatch,numcatch));%randomize order of catch trials within a block
trialAmountTotal = size(condMatrix,2);%new total trial amount

% Randomise the conditions
shuffler = Shuffle(1:trialAmountTotal);
condMatrixShuffled = condMatrix(:, shuffler);

% Make sure the same sequence never directly repeats (i.e. at least 1 other sequence in between)
for i=1:trialAmountTotal-2
    if condMatrixShuffled(1,i) == condMatrixShuffled(1,i+1)
        
        % if sequence i == sequence i+1, swap sequence i+1 with sequence i+2
        temp = condMatrixShuffled(:,i+1);
        condMatrixShuffled(:,i+1) = condMatrixShuffled(:,i+2);
        condMatrixShuffled(:,i+2) = temp;
    end
end

% add fourth row for type of catch trial (fixation vs. occlusion)
condMatrixShuffled(4,:) = NaN;

% Insert catch trials with a uniform distribution with the constraint:
% not too many normal trials in between two catch trials 
for icatch = 1:numcatch
    
    spot2insert = ceil(ceil((icatch-1)*((trialAmountTotal+numcatch)/numcatch)) + floor((trialAmountTotal+numcatch)/numcatch)*rand(1,1));
    
    condMatrixShuffled(:,spot2insert+1:trialAmountTotal+icatch) = condMatrixShuffled(:,spot2insert:trialAmountTotal+(icatch-1));
    condMatrixShuffled(1,spot2insert) = catchTrialPool{1,icatch};
    
    % condition:
    condMatrixShuffled(2,spot2insert) = randi(3,1);% any condition
    
    condMatrixShuffled(3,spot2insert) = icatch;
    
    condMatrixShuffled(4,spot2insert) = logical(catchTrialPool{3,icatch});% 0 = fixation, 1 = occlusion, NaN = not a catch trial
end
trialAmountTotal = trialAmountTotal + numcatch;

% MAKE AN OUTPUT MATRIX FOR STORING ALL VARIABLES PLUS RESPONSE DATA
% Matrix for storing behavioural data:
% 1st column is subject number
% 2nd column is block number
% 3nd column will record the trial number
% 4rd through 6th column trial information from CondMatrix
% 7th column correct versus incorrect for catch trials (i.e. 1 = correct, 0 = incorrect, 2 = not a catch trial)
% 8th column RT
% 9th column timing for comparing measured with intended timing
% 10th column measure actual ITI
% 11th column trigger value for that trial (in this experiment can simply be sequence number + catch trial info)
output = struct('subject_num',num2cell(repmat(subject_num,1,trialAmountTotal)),...
    'block_num',num2cell(repmat(block,1,trialAmountTotal)),...
    'trial_num',num2cell(1:trialAmountTotal),...
    'sequence',num2cell(condMatrixShuffled(1,:)),...
    'condition',num2cell(condMatrixShuffled(2,:)),...
    'catch',num2cell(condMatrixShuffled(3,:)),...
    'catchFixOrOcclude',num2cell(condMatrixShuffled(4,:)),...
    'correct_response',num2cell(NaN(1,trialAmountTotal)),...
    'RT',num2cell(NaN(1,trialAmountTotal)),...
    'testTiming',num2cell(NaN(1,trialAmountTotal)),...
    'testITI',num2cell(NaN(1,trialAmountTotal)),...
    'trigger',num2cell(NaN(1,trialAmountTotal)),...
    'time_catch',num2cell(NaN(1,trialAmountTotal))...
    );

% timings that are constant over trials
VidOnset = 0; % in seconds from trial start
VidOffset = 5;
maxTestTime = 3;% extra time to respond after video stops
timeCrossColor = 0.2;
timeOcclusion = 0.4;% time of occlusion

% timings that are random over trials
ITI = 1.8 + .4 * rand(1,trialAmountTotal);
for itrial = 1:trialAmountTotal
    output(itrial).ITI = ITI(itrial);
end

% create trigger values between 0 and 255 for MEG and eyetracker
% sequence onset trigger:
%       - 1:14 for sequence number
%       - + 100 for catch trial
%       - + 20 for condition 2, + 40 for condition 3
for itrial = 1:trialAmountTotal
    output(itrial).trigger = (output(itrial).condition-1)*20 + output(itrial).sequence + logical(output(itrial).catch) * 100;
end

%if not catch trial, response can already be set to 2
for itrial = 1:trialAmountTotal
    if ~output(itrial).catch
        output(itrial).correct_response = 2;
    end
end

%determine time of event in catch trials
catchcount = 0;
for itrial = 1:trialAmountTotal
    if ~output(itrial).catch
        output(itrial).time_catch = 0;
    else
        catchcount = catchcount + 1;
        output(itrial).time_catch = catchTrialPool{2,catchcount};
    end
end

if trackeye
    % INITIALIZE EYELINK
    % It is better not to send too many Eyelink commands to the eye-tracker in a row. For this reason, between them, we wait for a short time, here defined.
    elk.wait = 0.01;
    
    % This code initializes the connection with the eyelink: if something fails, it exit program with error
    if EyelinkInit()~= 1;
        error('Eyelink disconnected !!!');
    end;
    
    % We need to provide Eyelink with details about the graphics environment and perform some initializations. The initialization information is returned in a
    % structure that also contains useful defaults and control codes (e.g. tracker state bit and Eyelink key values). The structure, moreover, act asn an handle
    % for subsequent commands, like "windowHandle" for Psychtoolbox.
    elk.el = EyelinkInitDefaults(window);
    
    % Here we create the name for the eyelink datafile. Data gathered from the eye tracker are saved on the eye-tracking PC in a file. Data from all users are
    % saved in the same folder and the folder is routinely cleaned up without any advice. So, be sure to copy your data after the experiment and choose an
    % unique name for the datafile (containing date/time, subject number etc...). It has to be less than 8 characters long.
    %
    if subject_num < 10
        Eyefilename = ['ivUP0' num2str(subject_num) 'b' num2str(block)];
    else
        Eyefilename = ['ivUP' num2str(subject_num) 'b' num2str(block)];
    end
    Eyefilename = [Eyefilename '.edf'];
    elk.edfFile = sprintf(Eyefilename);		% Create file name
    Eyelink('Openfile', elk.edfFile);									% Open the file to the eye-tracker
    
    % Writing a short preamble to the file helps if the name became not that informative ;-)
    Eyelink('command', sprintf('add_file_preamble_text ''Ingmar de Vries Project Unpredict: subject %d ; block %d ; practice %d ; time %s''', subject_num, block, practice, datestr(now, 'YYYYmmddhhMM')));
    
    % Setting the eye-tracker so as to record GAZE of  LEFT and RIGHT eyes, together with pupil AREA
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
    
    % Setting the proper recording resolution, proper calibration type, as well as the data file content
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, screenXpixels - 1, screenYpixels - 1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, screenXpixels - 1, screenYpixels - 1);
    
    % Setting the proper calibration type. Usually we use 9 points calibration. For a long range mount also 13 points (HV13) is a good (longer) calibration.
    Eyelink('command', 'calibration_type = HV9');
    
    % Setting the proper data file content
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS');
    
    % Setting link data (used for gaze cursor, optional)
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');
    
    % Saccade detection thresholds (optional)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    
    % Now make sure that we are still connected to the Eyelink ... otherwise throw error
    if Eyelink('IsConnected')~=1
        error('Eyelink disconnected !!!');
    end;
    
    % EYELINK CALIBRATION
    % This code allow the EyeLink software to take control of your psychtoolbox screen. This means that at this point you will see participant eyes as recorded
    % by the Eye-tracker camera on the MEG whiteboard, a condition essential for setting up the camera. After setting up the camera you can perform calibration
    % and validation at this step.
    
    % Some calibration parameters
    elk.el.foregroundcolour = 0;
    elk.el.backgroundcolour = [bg bg bg] * 255;
    
    % Give eye-tracker control of the screen for camera setup and calibration, until you exit back to psychtoolbox by pressing ESC
    EyelinkDoTrackerSetup(elk.el);
    
end

if MEG
    % INITIALIZE DATAPIXX
    Datapixx('Open');					% Open DataPixx
    
    Datapixx('SetVideoMode', 0);		% This set video mode to normal passthrought, no stereo mode. C24, Straight passthrough from DVI 8-bit RGB to VGA RGB.
    % In this configuration luminance is linear with RGB (see our wiki).
    
    Datapixx('StopAllSchedules');		% Stop all schedules (audio waveforms, triggers etc...)
    
    Datapixx('SetDoutValues', 0);		% Set digital output to zero, as required to prepare for triggering
    
    Datapixx('EnableDinDebounce');		% Enable response debouncing. This is required to prune out spurious button presses after a real response
    
    Datapixx('SetDinLog');				% Clear digital input logger, i.e: clear old responses in the register
    Datapixx('StopDinLog');				% Stop running response logger
    
    Datapixx('RegWrRd');				% So far, no real changes occurred on the physical blue box devices. This command synchronize local and remote registers
    % in a read/write mode and immediately. Only now, the blue box status is as determined by the above initializations.
    
    responseButtonsMask = 2^0 + 2^1 + 2^2 + 2^3;	% Values of response buttons are stored in a cumbersome binary way. This is a binary mask useful to
    % transform them in decimal human-readable values. In particular, red = 1, blue = X, geen = X and yellow =
    % X. It works. Just believe it. I do, I am a true believer. Neo is the one.
end


% KEYBOARD (ONLY FOR BEHAVIOURAL VERSION)
KbName('UnifyKeyNames');
spaceKey = KbName('space');
escapeKey = KbName('ESCAPE');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
RestrictKeysForKbCheck([escapeKey spaceKey upKey downKey]);

%----------------------------------------------------------------------
%                       START WITH INSTRUCTIONS
%----------------------------------------------------------------------

% Instructions
if practice
    
    Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
    
    % Reset and fire up the response logger
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');                        % Commit changes to/from DP
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'Thank you for participating in this experiment',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to start the instructions -','center', yCenter+yCenter/2, black);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'Each trial starts with a fixation cross \n\n Whenever you see a fixation cross, please fixate it',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black);
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'Next you will see a video of a ballet dancer \n\n Pay attention to the dancer \n\n Keep fixating the cross while doing so',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black);
        Screen('DrawTexture', window, stimuli{1,1}, [], stimPos);
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'On some trials the fixation cross shortly becomes purple \n\n You have to indicate the color change \n\n by pressing the red button',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black);
        Screen('DrawTexture', window, stimuli{1,50}, [], stimPos);
        Screen('DrawLines', window, allCoords,lineWidthPix, colors(3,:), [xCenter yCenter], 2);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
        
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'If you need to respond \n\n you will receive direct feedback',...
            'center', yCenter/3, black);
        DrawFormattedText(window,...
            'Correct',...
            'center', yCenter*0.9, colors(1,:));
        DrawFormattedText(window,...
            'OR',...
            'center', yCenter, black);
        DrawFormattedText(window,...
            'Error',...
            'center', yCenter*1.1, colors(2,:));
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black)
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    DrawFormattedText(window,...
        'Please remember to relax and fixate the cross',...
        'center', 'center', black);
    DrawFormattedText(window,'- When ready, inform the experimenter -','center', yCenter+yCenter/2, black);
    Screen('Flip', window);
    KbStrokeWait;
    
else
    
    DrawFormattedText(window,...
        'I will now measure the position of your head \n\n Please do not move your head anymore from now on \n\n Please remember to relax and fixate the cross',...
        'center', 'center', black);
    DrawFormattedText(window,'- When ready, inform the experimenter -','center', yCenter+yCenter/2, black);
    Screen('Flip', window);
    KbStrokeWait;

end

%countdown at start each block
DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(3)], 'center', 'center', black);
Screen('Flip', window);
WaitSecs(1);
DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(2)], 'center', 'center', black);
Screen('Flip', window);
WaitSecs(1);
DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(1)], 'center', 'center', black);
Screen('Flip', window);
WaitSecs(1);

%----------------------------------------------------------------------
%                       START ACTUAL EXPERIMENT
%----------------------------------------------------------------------

% TRIAL LOOP
correctCount = 0;
catchcount = 0;
if practice
    trialAmountTotal = trialAmountPractice;
end

for trial = 1:trialAmountTotal
    
    % make sure there's a fixation cross in between trials while current
    % trial settings are being generated and eyelink recording is started
    Screen('FillRect', window, bg, [xCenter-(screenXpixels/2) yCenter-(screenYpixels/2) xCenter+(screenXpixels/2) yCenter+(screenYpixels/2)]);
    Screen('FillRect', window, colorPlaceHolder, [xCenter-(stimWidth/2) yCenter-(stimHeight/2) xCenter+(stimWidth/2) yCenter+(stimHeight/2)]);
    Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
    
    %photo diode
    Screen('FillRect', window, [0 0 0], [0 0 30 30]);
    
    ITIstart = Screen('Flip', window);
    
    % initialize some things and extract from output struct (saves time to do it here and not while video is being shown)
    response=0;
    respTime = 0;
    rtMarkerTime = 0;
%     testTime = 0;
    currentSeq = output(trial).sequence;
    currentCond = output(trial).condition;
    currentCatch = output(trial).catch;
    currentStim = squeeze(stimuli(currentCond,currentSeq,:));
    currentCatchTime = output(trial).time_catch;
    
    if currentCatch
        postOcclusionStim = catchTrialPool{3,currentCatch};% 0 = fixation cross catch trial
        currentCorrectResp = (2 - catchTrialPool{4,currentCatch})^3;%1 = correct = red, 8 = incorrect = blue
        if postOcclusionStim ~= 0% if not fixation cross trial
            if currentCond == 3% in scrambled condition, change stimulus to 'normal' from 500 msec before occlusion so catch trial is not impossibly hard
                currentStim(round((currentCatchTime-.5)*refrate):end) = squeeze(stimuli(1,currentSeq,round((currentCatchTime-.5)*refrate):end));
                
                %change to different trial after occlusion if it's supposed to be a mismatch, or keep same if it's supposed to be a 'match' catch trial
                currentStim(round((currentCatchTime+timeOcclusion)*refrate):end) = squeeze(stimuli(1,postOcclusionStim,round((currentCatchTime+timeOcclusion)*refrate):end));
            else
                %change to different trial after occlusion if it's supposed to be a mismatch, or keep same if it's supposed to be a 'match' catch trial
                currentStim(round((currentCatchTime+timeOcclusion)*refrate):end) = squeeze(stimuli(currentCond,postOcclusionStim,round((currentCatchTime+timeOcclusion)*refrate):end));                
            end

        end
    
    end
        
    freezeFrame = currentCatchTime * refrate;
    
    if trackeye
        % EYELINK RECORDING
        
        Eyelink('Message', 'TRIALID %d', trial);
        Eyelink('command', 'record_status_message "SUBJECT = %d ; BLOCK = %d ; TRIAL = %d"', subject_num, block, trial);
        
        % As specified before, it is better not to send to many Eyelink commands to the eye-tracker in a row.  So, after the command, we wait the pre-set time
        WaitSecs(elk.wait);
        
        % Here we start recording eyelink data (left/right gaze and pupil size), preceded by a short pause
        Eyelink('Command', 'set_idle_mode');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %transfer image to host
        [width, height]=Screen('WindowSize', screenNumber);
        imgfile= ([rootdir '\experiment\fixation.bmp']);
        transferimginfo=imfinfo(imgfile);
        
        fprintf('img file name is %s\n',transferimginfo.Filename);
        
        % image file should be 24bit or 32bit bitmap
        % parameters of ImageTransfer:
        % imagePath, xPosition, yPosition, width, height, trackerXPosition, trackerYPosition, xferoptions
        transferStatus =  Eyelink('ImageTransfer',transferimginfo.Filename,0,0,transferimginfo.Width,transferimginfo.Height,width/2-transferimginfo.Width/2 ,height/2-transferimginfo.Height/2,1);
        if transferStatus ~= 0
            fprintf('*****Image transfer Failed*****-------\n');
        end
        
        WaitSecs(0.1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%
        WaitSecs(elk.wait);
        Eyelink('StartRecording', 1, 1, 1, 1);
        %%
        WaitSecs(elk.wait);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    % Continue with fixation cross in between trials (ITI)
    Screen('FillRect', window, bg, [xCenter-(screenXpixels/2) yCenter-(screenYpixels/2) xCenter+(screenXpixels/2) yCenter+(screenYpixels/2)]);
    Screen('FillRect', window, colorPlaceHolder, [xCenter-(stimWidth/2) yCenter-(stimHeight/2) xCenter+(stimWidth/2) yCenter+(stimHeight/2)]);
    Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
    %photo diode
    Screen('FillRect', window, [0 0 0], [0 0 30 30]);
    
    ITImiddle = Screen('Flip', window);
    WaitSecs(output(trial).ITI - (ITImiddle - ITIstart));
    % Note that use of WaitSecs is not very accurate,
    % but for the ITI it does not matter that much. It's randomly jittered anyway.
    % What is crucial is the timing of sequence onset and each frame in the sequence
    % That is why we do this with counting exact frames (see below)
    testITI2 = GetSecs;
    output(trial).testITI = testITI2-ITIstart;
    %------------------------------------------------------------------
    %              START ACTUAL TRIAL WITH ACCURATE TIMING
    %------------------------------------------------------------------
    
    disco = 0;
    for frame = 1:ceil((VidOffset+maxTestTime) * refrate)
        
        %always draw black placeholder at location of video
        Screen('FillRect', window, colorPlaceHolder, [xCenter-(stimWidth/2) yCenter-(stimHeight/2) xCenter+(stimWidth/2) yCenter+(stimHeight/2)]);
        
        %duration of video
        if frame <= ceil(VidOffset * refrate)
            
            %photo diode
            disco = 1 - disco;%during duration of video let square for photodiode flicker between black and white
            Screen('FillRect', window, [disco disco disco], [0 0 30 30]);
            
            % Draw each frame of the video
            Screen('DrawTexture', window, currentStim{frame}, [], stimPos);
        end
        
        % Fixation cross
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        
        % if catch trial, and currentCatchType = 0, fixation cross shortly turns purple
        if currentCatch && postOcclusionStim == 0 && (frame > ceil(currentCatchTime * refrate)) && (frame <= ceil((currentCatchTime+timeCrossColor) * refrate))
            Screen('DrawLines', window, allCoords,lineWidthPix, colors(3,:), [xCenter yCenter], 2);
            
            % if catch trial, and currentCatchType is not 0, then there's the occlusion
        elseif currentCatch && postOcclusionStim~=0 && (frame > ceil(currentCatchTime * refrate)) && (frame <= ceil((currentCatchTime+timeOcclusion) * refrate))
            % draw occlusion at location of video
            Screen('FillRect', window, white-.2, [xCenter-(stimWidth/2) yCenter-(stimHeight/2) xCenter+(stimWidth/2) yCenter+(stimHeight/2)]);
            
        end
        
        % TRIGGERS
        % memory display onset
        if MEG
            if frame == 1
                
                %Datapixx triggering
                triggerPulse = [1 0] .* (output(trial).trigger);
                Datapixx('StopDoutSchedule');
                Datapixx('WriteDoutBuffer', triggerPulse);
                Datapixx('SetDoutSchedule', 1.0/refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                Datapixx('StartDoutSchedule');
                
                % Eyelink triggering
                Eyelink('Message', sprintf('Video onset %d', output(trial).trigger));
                
                % White square for photodiode
                Screen('FillRect', window, [1 1 1], [0 0 30 30]);
                            
                Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                
            elseif currentCatch && frame == ceil(currentCatchTime * refrate)+1
                                
                %Datapixx triggering
                triggerPulse = [1 0] .* (output(trial).trigger+100);
                Datapixx('StopDoutSchedule');
                Datapixx('WriteDoutBuffer', triggerPulse);
                Datapixx('SetDoutSchedule', 1.0/refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                Datapixx('StartDoutSchedule');
                
                % Eyelink triggering
                Eyelink('Message', sprintf('Video onset %d', output(trial).trigger+100));
                
                Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                
                % Reset and fire up the response logger
                Datapixx('SetDinLog');
                Datapixx('StartDinLog');
                
            end
            
            %Now send all instructions to the Datapixx box, as close as possible to the actual screen flip in PTB
            Datapixx('RegWrVideoSync');
            
        end
        
        Screen('Flip',window);
        
        if currentCatch && frame == ceil(currentCatchTime * refrate)+1
            % This has to be in its own loop after videosync and
            % screenflip, otherwise causes a delay
            Datapixx('SetMarker');
            Datapixx('RegWrRd');
            rtMarkerTime = Datapixx('GetMarker');
        end
        
        % in catch trials a response is required
        if currentCatch && (frame > ceil(currentCatchTime * refrate))
            
            if MEG
                Datapixx('RegWrRd');						% Commit changes to/from DP
                status = Datapixx('GetDinStatus');			% Get response logger status
                
                if status.newLogFrames > 0					% We've got new data in response buffer !!!
                    
                    [data, time] = Datapixx('ReadDinLog');	% Read data in
                    
                    response = bitand(data(end), responseButtonsMask);
                    respTime = time(end);	%
                    
                    % Eyelink triggering
                    Eyelink('Message', sprintf('Response %d', response));
                    
                    % Send RT/response trigger (response value + 128)
                    Datapixx('EnableDinDebounce');
                    Datapixx('StopDoutSchedule');
                    triggerPulse = [1 0] .* (response+60);
                    Datapixx('WriteDoutBuffer', triggerPulse);
                    Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                    Datapixx('StartDoutSchedule');
                    Datapixx('RegWr');
                    
                    break
                    
                end
                
            else
                [keyIsDown,respTime,keyCode] = KbCheck;
                if keyCode(upKey)
                    response = 1;
                    if trackeye
                        Eyelink('Message', sprintf('Response %d', response));
                    end
                    break
                elseif keyCode(downKey)
                    response = 2;
                    if trackeye
                        Eyelink('Message', sprintf('Response %d', response));
                    end
                    break
                end
            end
            
        end
        
        %temporary timing test
        if frame == 1
            test1 = GetSecs;
        end
        if frame == ceil(VidOffset*refrate)
            test2 = GetSecs;
            output(trial).testTiming = test2-test1;
        end
        
        if currentCatch && (frame == ceil(currentCatchTime * refrate))
            testTime = GetSecs;
        end
        
        % check for esc key
        [keyIsDown,respTime,keyCode] = KbCheck;
        if keyCode(escapeKey)
            sca;
            disp('The escape key has been pressed');
            Datapixx('Close');
            return
        end
        
        if ~currentCatch && (frame > ceil(VidOffset * refrate))
            break
        end
    end
    
    if trackeye
        % Stop eyelink
        Eyelink('StopRecording');
    end
    
    if currentCatch
        % check RT
        output(trial).RT = respTime - rtMarkerTime;
        % respTime and testTime set to 0 at start each trial
        % if respTime is still 0 here no response was given on this trial
        % if testTime is still 0 here this trial was not a catch trial
        
        % feedback
        Screen('FillRect', window, bg, [xCenter-(screenXpixels/2) yCenter-(screenYpixels/2) xCenter+(screenXpixels/2) yCenter+(screenYpixels/2)]);
        
        if response == currentCorrectResp
            
            output(trial).correct_response = 1;
            correctCount = correctCount + 1;
            
            DrawFormattedText(window, 'Correct', 'center', 'center', colors(1,:));
            
        elseif (9-response == currentCorrectResp) && postOcclusionStim~=0
            
            output(trial).correct_response = 0;
            if practice
                DrawFormattedText(window,...
                    'Error \n\n Pay attention to the dancer \n\n while fixating the cross',...
                    'center', 'center', colors(2,:));
            else
                DrawFormattedText(window,...
                    'Error',...
                    'center', 'center', colors(2,:));
            end
                        
        elseif response == 0 && postOcclusionStim==0
            
            output(trial).correct_response = 0;
            DrawFormattedText(window,...
                'Missed the color change in the fixation cross!',...
                'center', 'center', colors(2,:));
            
        elseif response == 0 && postOcclusionStim~=0
            
            output(trial).correct_response = 0;
            DrawFormattedText(window,...
                'Too slow, please respond faster next time',...
                'center', 'center', colors(2,:));
            
        end
        
        disp(['response :' num2str(response) ' , and RT: ' num2str(output(trial).RT)]);
                
        Screen('Flip', window);
        WaitSecs(2);
    end
    
    %Save each trial so that data is stored if the experiment stops
    %unexpectedly
    save([savedir filename], 'output');
    
end

% FEEDBACK AFTER EACH BLOCK
percentCorrect = ceil(correctCount/(sum(logical(condMatrixShuffled(3,1:trialAmountTotal))))*100);

DrawFormattedText(window,['You had ' num2str(percentCorrect) '% correct in this block'],'center', yCenter/2, black);
DrawFormattedText(window,'Please tell the experimenter that you finished this block \n\n and take a small break','center', yCenter+yCenter/2, black);
Screen('Flip', window);
KbStrokeWait;

% CLOSE STUFF
% Eyelink
% Download data file
if trackeye
    try
        fprintf('Receiving data file ''%s''\n', elk.edfFile );
        status = Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2 == exist(elk.edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', elk.edfFile, savedir );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', elk.edfFile );
    end
    
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    Eyelink('ShutDown');
end

if MEG
    % Close DataPixx
    Datapixx('StopAllSchedules');		% Stop all schedules
    Datapixx('SetDoutValues', 0);       % Reset triggers to zero
    Datapixx('StopDinLog');             % Stop response buttons recording
    Datapixx('Close');                  % Close DataPixx
end

% Psychtoolbox
Priority(0);
sca

disp(['You had ' num2str(percentCorrect) '% correct in this block']);
if practice
    disp(['Refresh rate is ' num2str(refrate)]);
end

if block == 1
    save(fn_backup,'samples_all','randID_all');
end