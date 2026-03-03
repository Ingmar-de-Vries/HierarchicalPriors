function [vecrepGRAYSMOOTH] = HierarchicalPriors_video2vector(cfg)
%% Transform video frames to a low level vector representation for each frame
%INPUT
% video = full name (incl path) of video

videoHeader = VideoReader(cfg.videoName);
video = read(videoHeader);

%some preprocessing of videos:
video = im2double(video);% correctly transform uint8 to double

% Vector representation of CIELAB values for each frame
videoLAB = rgb2lab(video);
clear video

% Vector representation of grayscale / luminance values approximated by the L (first) dimension in the CIELAB colorspace 
videoGRAY = squeeze(videoLAB(:,:,1,:));
clear videoLAB

% Vector representation of grayscale smoothed with Gaussian kernel
% Kriegeskorte 2008 FHN uses Gaussian kernel with FWHM = 11.75 pixels
% Conversion to sigma is roughly: sigma = FWHM / 2.355, i.e.: 11.75/2.355 ~= 5
videoGRAYSMOOTH = imgaussfilt(videoGRAY,5);
clear videoGRAY

% reduce file size and reshape
% videoGRAYSMOOTH = uint8(rescale(videoGRAYSMOOTH,1,2^8));
vecrepGRAYSMOOTH = reshape(videoGRAYSMOOTH,numel(videoGRAYSMOOTH(:,:,1)),size(videoGRAYSMOOTH,3));

end