%% diffuserRecon_fig4v1.m
% Neerja Aggarwal & Joseph Malone
% July 13th, 2022
% Purpose: generate the recon for measurements.  

%% load calibration dataset (takes a minute)
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
spectralPSF_3D = tifRead(filename);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description

%% TODO: % make a function that creates a sampling mask. takes in three parameters:
% example imagesize, min square,max square, number of sampling points. 
% use randSamplePSF then ignore any values inside min square.

maskbox = 500; %half width of the box. 

samplepercent = 0.1;

[N1,N2,N3] = size(spectralPSF_3D);
sampfac = round(sqrt(N1*N2*(samplepercent/100)));

samp_xy(:,1) = randi(N1,[sampfac,1]);
samp_xy(:,2) = randi(N2,[sampfac,1]);

% compute portion inside the box
logic(:,1) = samp_xy(:,1)>(floor(N1/2)-maskbox);
logic(:,2) = samp_xy(:,1)<(floor(N1/2)+maskbox);
logic(:,3) = samp_xy(:,2)>(floor(N2/2)-maskbox);
logic(:,4) = samp_xy(:,2)<(floor(N2/2)+maskbox);
% remove the indices that are inside the box
insidebox = all(logic,2);

%% use linear indexing instead

% 1. get random samples over linear indexing. 
maskbox = 500; %half width of the box. 
samplepercent = 0.1;
[N1,N2,N3] = size(spectralPSF_3D);
numsamples = floor(N1*N2*(samplepercent/100));
samp_ind = randi(N1*N2,[numsamples,1]);

% 2. get indices for mask
rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];
maskind = sub2ind([N1 N2],rows,cols)';

%% set value to NaN inside A for mask
spectralPSF_3D(rows,cols,:) = nan;

s = spectralPSF_3D(:,:,1);

for m = 1:N3
    %reshape 
    s = spectralPSF_3D(:,:,m);
    samp_full(:,m) = s(samp_ind); %use linear indexing for truly random

nanind = isnan(samp_full);
% remove those values here and in samp_ind
spectralPSF_2D = samp_full(~nanind);

% 4. get samples of A. 



%%
samp_full = spectralPSF_3D(samp_xy(:,1),samp_xy(:,2),:);
spectralPSF_2D = reshape(samp_full,[sampfac^2,N3]);



% Option 2: randomly sampling.  Then throw away the values in the mask
% range. compute the portion thrown away.  resample and repeat until the
% samp_xy list is long enough? 

