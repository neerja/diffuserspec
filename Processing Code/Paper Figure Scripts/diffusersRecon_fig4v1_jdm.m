%% Figure 4 processing - broadband diffuser recon comparison

% SET DIR to - Projects\smartOCT\Analysis\Data\CDMRPVRP\diffuserSpecData_OpticsLetters\organizedCodeData\Data
disp('select main data folder')
originalDir = uigetdir;

%% load in pre calibrated data and constants 

% select calibrationFiles folder (located with
disp('select calibration folder')
calibrationDir = uigetdir;
cd(calibrationDir)

load('calibrationOrig.mat');
load('calibrationCurve.mat');
load('signalRange.mat');
load('calibrationWavelengths_fit.mat');
load('DOE_broadband_source1.mat');
load('DOE_wavelengths.mat');
load('tape_broadband_source1.mat')
load('tape_wavelengths.mat');
load('samp_ind_nonan.mat');
cd(originalDir)

%% load diffuser spec data

%load calibration matrix
disp('Loading calibration matrix')
SSTM = uigetfile('*tif');
spectralPSF_3D = tifRead(SSTM);

%load background frame
disp('Loading background signal for background subtraction')
SSTM_bg = uigetfile('*tif');
backgroundFrame = tifRead(SSTM_bg);
%subtract background and resign negative values
spectralPSF_3D = abs(spectralPSF_3D - backgroundFrame);

% get average intensity for each spectral PSF and normalize max to 1
normalizedIntensity = mean(mean(spectralPSF_3D,1),2);
normalizedIntensity = normalizedIntensity/max(normalizedIntensity);

% normalize spectral intensity and set max value to 1
spectralPSF_3D = spectralPSF_3D./normalizedIntensity;
spectralPSF_3D = spectralPSF_3D/max(spectralPSF_3D(:));

clear backgroundFrame normalizedIntensity

%% load broadband data

% select calibrationFiles folder (located with
calibrationDir = uigetdir;
cd(calibrationDir)

disp('select broadband diffuse spectrum')
broadband_fname = uigetfile('*tif');
broadband = mean(tifRead(broadband_fname),3);

disp('select broadband diffuse spectrum background frame')
broadbandbg_fname = uigetfile('*tif');
broadband_bg = mean(tifRead(broadbandbg_fname),3);

broadband = abs(broadband-broadband_bg);

clear broadband_bg broadband_fname broadbandbg_fname

cd(originalDir)

%% get sampling parameters

samplePercent = .1;
maskSize = 200;
% [spectralPSF_2D, samp_ind_nonan] = ...
%     randSamplePSF_mask_jdm(spectralPSF_3D,samplePercent,maskSize);
spectralPSF_3D = squeeze(reshape(spectralPSF_3D,[],1,size(spectralPSF_3D,3)));
spectralPSF_2D = spectralPSF_3D(samp_ind_nonan,:);

broadband_reshape = squeeze(reshape(broadband,[],1));
broadband_forRecon = broadband_reshape(samp_ind_nonan,:);

%%
sig =12;
[recon_broadband_gSVD] = gaussSVD(spectralPSF_2D(:,1:344),...
    broadband_forRecon,sig);
recon_broadband_gSVD = abs(recon_broadband_gSVD/max(recon_broadband_gSVD(:)));

figure; plot(calibrationWavelengths_fit,recon_broadband_gSVD);hold on
plot(full_gt_waves(32:546),full_gt(32:546)/max(full_gt(32:546)));
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')
legend('Recon','Ground Truth')