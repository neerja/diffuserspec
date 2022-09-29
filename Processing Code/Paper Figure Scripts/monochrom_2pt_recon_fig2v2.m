% monochrom_2pt_recon_fig2v2.m
% Neerja Aggarwal
% July 20th, 2022

% Use gauss SVD to make monochromatic recon. Uses random sampling with box
% in middle. 

%% load calibration dataset (takes a minute)
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
spectralPSF_3D = tifRead(filename);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description
[N1,N2,N3] = size(spectralPSF_3D);
% Create randomly sampled 2D spectralPSF
% choose percentage of data to sample
samplepercent = 0.1; %percent of pixels
maskbox = 500; %radius to ignore in center
[spectralPSFrand, samp_ind_nonan] = randSamplePSF_mask(spectralPSF_3D,samplepercent,maskbox);

%% Figure 1a - monochromatic reconstruction

% down sample 2D transfer function by a factor of 2 for reconstruction, code 
% takes every other frame starting at 1 (all odd frames)
    
downSampFac = 2;
spectralPSF_2D_downSamp= spectralPSFrand(:,1:downSampFac:end);
numFrames = 6;
wv_downSamp = calibrationWavelengths_fit(1:downSampFac:end);    
% grab 6 monochromatic frames evenly distributed over full bandwidth
% note: each frame is an even indexed frame and not part of the downsampled
% transfer matrix

stepsize = ceil(N3/numFrames);
offset = floor(round(N3/numFrames)-(round(N3/numFrames)/2));

multiPeakVec = [stepsize-offset,(stepsize*2)-offset,(stepsize*3)-offset,...
    (stepsize*4)-offset,(stepsize*5)-offset,(stepsize*6)-offset];
spectrumForRecon_multiWave = spectralPSFrand(:,multiPeakVec);
    
% set threshold for Truncated SVD recon and run recon and set max value of
% recon to 1
thresh = 16;
[recon] = TruncSVD(spectralPSF_2D_downSamp,...
        spectrumForRecon_multiWave,thresh);
 recon_full = abs(recon/max(recon(:)));

figure;
plot(wv_downSamp',recon_full)

%% try with Gauss recon

sigma = 20;
[recon] = gaussSVD(spectralPSF_2D_downSamp,spectrumForRecon_multiWave,sigma);
 recon_full = abs(recon/max(recon(:)));

figure;
plot(wv_downSamp',recon_full)
xlabel('Wavelength (nm)')
ylabel('Intensity')

% list ground truth wavelengths and compare to peak value:

gt = calibrationWavelengths_fit(multiPeakVec)
[m,ind] = max(recon_full);
mp = wv_downSamp(ind)

error = mean(abs(gt - mp))

% accuracy is +- 0.25 nm. 

%% FIGURE 1b two-point resolution
    
% select frames for 1nm and 2nm spacing (each wavelength step,
% corresponds to .25 nm so 1nm step = 4 frame seperation, 2nm step = 8
% frame seperation. Frames are averaged together to simulate
% simultaneous acquisition

nm2_ind = [172 180];
nm1p5_ind = [172 188];

res2nm = mean(spectralPSFrand(:,nm2_ind),2);
res15nm = mean(spectralPSFrand(:,nm1p5_ind),2);
    
sigma = 23;
[recon] = gaussSVD(spectralPSF_2D_downSamp,res15nm,sigma);
 recon_full = abs(recon/max(recon(:)));

figure;
plot(wv_downSamp',recon_full)
xlabel('Wavelength (nm)')
ylabel('Intensity')