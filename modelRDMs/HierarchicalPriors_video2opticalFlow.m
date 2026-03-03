function [vecrepDIR, vecrepMAG] = HierarchicalPriors_video2opticalFlow(cfg)
%% Transform video frames to an optical flow vector representation for each frame
%INPUT
% video = full name (incl path) of video

% estimate optical flow of each frame
opticFlow = opticalFlowFarneback;

% load video
vidReader = VideoReader(cfg.videoName);

% initialize variables
matrepMAG = zeros(vidReader.NumFrames,vidReader.Height,vidReader.Width);
matrepORI = zeros(vidReader.NumFrames,vidReader.Height,vidReader.Width);
clear flowFB
iframe = 0;
while hasFrame(vidReader)
    
    iframe = iframe+1;
    
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    
    % estimate optical flow
    flowFB = estimateFlow(opticFlow,frameGray);
    
    matrepMAG(iframe,:,:) = flowFB.Magnitude;
    matrepORI(iframe,:,:) = flowFB.Orientation;
    
end

%tNew and t4deriv are only used for interpolation to correct time stamps
sampfreq = 50;
tNew = 0:1/sampfreq:5-1/sampfreq;
t4deriv = tNew+(1/sampfreq)/2;% timesteps in between for any derivative (e.g., position --> motion, or motion --> acceleration)
t4deriv(end) = [];

% interpolate because optical flow vectors now placed in between time points
matrepMAG = interp1(t4deriv,matrepMAG(2:end,:,:),tNew,'pchip','extrap');

% we need to unwrap the angles before we can interpolate them
anguw=unwrap(matrepORI);
anginterp=interp1(t4deriv,anguw(2:end,:,:),tNew,'pchip','extrap');
anginterp=mod(anginterp,2*pi);
anginterp(anginterp>pi)=anginterp(anginterp>pi)-2*pi;
matrepDIR(:,1,:,:) = cos(anginterp);
matrepDIR(:,2,:,:) = sin(anginterp);

%% last steps before creating vector, following Efros et al. 2003 
% spatially smooth with Gaussian filter
sigma = 5;% ~ as Kriegeskorte 2008 for pixelwise similarity
for iframe = 1:size(matrepDIR,1)
    
    % magnitude
    matrepMAG(iframe,:,:) = imgaussfilt(squeeze(matrepMAG(iframe,:,:)),sigma);
    
end

% store optical flow vector representation of each video
vecrepDIR = reshape(matrepDIR,size(matrepDIR,1),numel(squeeze(matrepDIR(1,:,:,:))));
vecrepMAG = reshape(matrepMAG,size(matrepMAG,1),numel(squeeze(matrepMAG(1,:,:,:))));

end