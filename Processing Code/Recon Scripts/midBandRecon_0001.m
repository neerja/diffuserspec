% midBandRecon_0001.m
% Neerja Aggarwal
% June 28th, 2022

% Purpose: try different algorithms for midband data recon on _#0001.  Clean up Joe's
% processing code to rely on reusable helper functions. 

% MAKE SURE running from repo base directory: diffuserspec/

%% load calibration dataset (takes a minute)
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
spectralPSF_3D = tifRead(filename);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description
%% Create randomly sampled 2D spectralPSF
% choose percentage of data to sample
samplepercent = 0.1; %percent of pixels
[spectralPSF_2D,sampxy,sampfac] = randSamplePSF(spectralPSF_3D,samplepercent);

%% load the test measurement
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0001.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';

spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);
%% load gt spectrum 
% load calibration files  WHAT ARE THESE?
load './Datasets matFiles/120grit/wavelength_gt.mat'

filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum_gt/partialSpectrum_1.spf2';

[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

figure;
hold on
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt)
xlabel('Wavelength (nm)')
ylabel('normallized intensity')
legend('groundtruth')

%% Try truncated SVD
thresh = 0.01; %2 is optimal
[recon,Alr,Ainv,Dprime,Dprime_thresh] = TruncSVD(spectralPSF_2D,...
    spectrumForRecon_sampled,thresh);
recon = abs(recon/max(recon(:)));
figure;
plot(calibrationWavelengths_fit,recon)
hold on
plot(calibrationWavelengths_fit,partialSpectrum_gt,'--')

xlabel('Wavelength (nm)')
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend('Reconstructed','Groundtruth')

%% Try truncated SVD  with the gaussian filter
sigma = 30;
recon_gauss = gaussSVD(spectralPSF_2D,spectrumForRecon_sampled,sigma);
fig = figure(1);
plot(calibrationWavelengths_fit,recon_gauss)
hold on
% plot(calibrationWavelengths_fit,spectrum_gt,'--')
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt)

xlabel('Wavelength (nm)')
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend('Reconstructed','Groundtruth')

% save figure in results
% 
% resultfile = './Results/120grit_midband0001_gaussSVD_30';
% saveas(fig, resultfile,'png')
% saveas(fig, resultfile,'eps')

%% Try ridge regression
lambda = 10000;
% matlab statistical and machine learning toolbox has builtin ridge
% regression
recon_ridge = ridge(spectrumForRecon_sampled, spectralPSF_2D,lambda);
recon_ridge = recon_ridge/max(recon_ridge);
fig = figure(1);
plot(calibrationWavelengths_fit,recon_ridge)
hold on
% plot(calibrationWavelengths_fit,spectrum_gt,'--')
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt)

xlabel('Wavelength (nm)')
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend('Reconstructed','Groundtruth')

