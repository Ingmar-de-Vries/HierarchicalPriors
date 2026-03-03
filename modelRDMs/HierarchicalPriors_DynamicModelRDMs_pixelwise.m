function HierarchicalPriors_DynamicModelRDMs_pixelwise(cfg)

% create dynamic RDMs from pixelwise grayscale of videos

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

% Loop over each combination of videos to compute RDMs
frames = 250;
newframes = 500;
RDMgraysmooth = zeros(length(vidNames),length(vidNames),frames,frames);
for iMovie1=1:length(vidNames)
     
    %Call function to extract low level vector representations
    IdV_cfg = [];
    IdV_cfg.videoName = fullfile(StimDir,vidNames(iMovie1).name);
    vecrepGRAYSMOOTH1 = HierarchicalPriors_video2vector(IdV_cfg);
    
    for iMovie2=1:length(vidNames)
        
        IdV_cfg = [];
        IdV_cfg.videoName = fullfile(StimDir,vidNames(iMovie2).name);
        vecrepGRAYSMOOTH2 = HierarchicalPriors_video2vector(IdV_cfg);

        clc;
        disp(['Correlating video ' num2str(iMovie1) ' with ' num2str(iMovie2)]);    
        
        RDMgraysmooth(iMovie1,iMovie2,:,:) = 1 - corr(vecrepGRAYSMOOTH1,vecrepGRAYSMOOTH2);
        
    end
end

% interpolate to 100Hz to match neural and kinematic RDMs in main script
fsVid = 50;
fsNew = 100;
tVid = 0:1/fsVid:5-1/fsVid;
tNew = 0:1/fsNew:5-1/fsNew;
tempall = zeros(length(vidNames),length(vidNames),newframes,newframes);
for istim1 = 1:length(vidNames)
    for istim2 = 1:length(vidNames)
        
        temp = squeeze(RDMgraysmooth(istim1,istim2,:,:));
        % bilinear interpolation:
        temp = interp1(tVid,temp,tNew,'pchip','extrap');
        temp = interp1(tVid,temp',tNew,'pchip','extrap')';
        tempall(istim1,istim2,:,:) = temp;
        
    end% second stimuli loop
end% first stimuli loop

RDMgraysmooth = tempall;

% save correlation matrices
save([DataDir filesep 'HierarchicalPriorsn_dynRDM_graysmooth'],'RDMgraysmooth');

end