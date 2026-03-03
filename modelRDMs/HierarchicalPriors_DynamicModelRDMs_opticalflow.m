function HierarchicalPriors_DynamicModelRDMs_opticalflow(cfg,~,~,~)

% create dynamic RDMs from optical flow model of videos

% set directories
if isfield(cfg,'cluster')
    addpath('//XXX/HierarchicalPriors/toolboxes/fieldtrip-20191113');
    addpath('//XXX/HierarchicalPriors/code/modelRDMs');
    rootdir = '//XXX/HierarchicalPriors';
else
    addpath('\\XXX\HierarchicalPriors\toolboxes\fieldtrip-20191113');
    addpath('\\XXX\HierarchicalPriors\code\modelRDMs');
    rootdir = '\\XXX\HierarchicalPriors';
end

% path to the movie data (mp4 or mov)
StimDir = fullfile(rootdir,'experiment','videos');
DataDir = fullfile(rootdir,'data','modelRDMs');

vidNames = dir(fullfile(StimDir, '*vid.mp4'));

%% first load videos, create vector representation of optical flow field, and store vector representation to disk
computeOFV = 0;
if computeOFV == 1
    for imov = 1:length(vidNames)

        clc;
        disp(['Creating vector representation of optical flow field for video ' num2str(imov)]);

        M1 = [];
        cfg_in.videoName = fullfile(StimDir,vidNames(imov).name);
        [M1.vecrepdir, M1.vecrepmag] = HierarchicalPriors_video2opticalFlow(cfg_in);

        save(fullfile(StimDir,['HierarchicalPriors_OF_vid' vidNames(imov).name(end-8:end-7)]),'M1','-v7.3');

    end
end

%% Next load vector representations and compute RDMs 
frames = 250;
newframes = 500;
RDMoptflow_dir = zeros(length(vidNames),length(vidNames),frames,frames);
RDMoptflow_mag = zeros(length(vidNames),length(vidNames),frames,frames);
for imov1 = 1:length(vidNames)
    
    % load first video
    load(fullfile(StimDir,['HierarchicalPriors_OF_vid' vidNames(imov1).name(end-8:end-7)]),'M1');
    mov1 = M1;
    
    for imov2 = 1:length(vidNames)
        
        clc;
        disp(['Correlation optical flow representation video ' num2str(imov1) ' with video ' num2str(imov2)]);
        
        % load second video
        load(fullfile(StimDir,['HierarchicalPriors_OF_vid' vidNames(imov2).name(end-8:end-7)]),'M1');
        M2 = M1;
        clear M1
        
        RDMoptflow_dir(imov1,imov2,:,:) = 1 - corr(mov1.vecrepdir',M2.vecrepdir');
        
        RDMoptflow_mag(imov1,imov2,:,:) = 1 - corr(mov1.vecrepmag',M2.vecrepmag');
        
    end
    
end

% interpolate to 100Hz to match neural and kinematic RDMs in main script
fsVid = 50;
fsNew = 100;
tVid = 0:1/fsVid:5-1/fsVid;
tNew = 0:1/fsNew:5-1/fsNew;
tempall1 = zeros(length(vidNames),length(vidNames),newframes,newframes);
tempall2 = tempall1;
for istim1 = 1:length(vidNames)
    for istim2 = 1:length(vidNames)
        
        temp = squeeze(RDMoptflow_dir(istim1,istim2,:,:));
        % bilinear interpolation:
        temp = interp1(tVid,temp,tNew,'pchip','extrap');
        temp = interp1(tVid,temp',tNew,'pchip','extrap')';
        tempall1(istim1,istim2,:,:) = temp;
        
        temp = squeeze(RDMoptflow_mag(istim1,istim2,:,:));
        % bilinear interpolation:
        temp = interp1(tVid,temp,tNew,'pchip','extrap');
        temp = interp1(tVid,temp',tNew,'pchip','extrap')';
        tempall2(istim1,istim2,:,:) = temp;
        
    end% second stimuli loop
end% first stimuli loop

RDMoptflow_dir = tempall1;
RDMoptflow_mag = tempall2;


save([DataDir filesep 'HierarchicalPriors_dynRDM_OFdir'],'RDMoptflow_dir');
save([DataDir filesep 'HierarchicalPriors_dynRDM_OFmag'],'RDMoptflow_mag');

