%% show measurement
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0002.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';

filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0002.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';

spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename);

a = min(min(spectrumForRecon_sampled));
b = max(max(spectrumForRecon_sampled));
bscale = .5;
figure;
imagesc(spectrumForRecon_sampled,[a,b*bscale]); %contrast stretched
box off;
axis off;
colormap('gray')

%% show spectral3Dtif 

filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
spectralPSF_3D = tifRead(filename);

% choose n frames

framenum =  [1,2,343,344] %choose which frames to plot

stm = spectralPSF_3D(:,:,framenum);

clear spectralPSF_3D
%% make frames for calibration
close all;

a = min(min(min(stm)));
b = max(max(max(stm)));
bscale = 0.1;
clim = [a, b*bscale];


for k1 = 1:length(framenum)
    y = squeeze(stm(:,:,k1));
    figure;
    imagesc(y,clim); %contrast stretched
    % colormap([1,0,0]);
    colormap('gray');
    box off;
    axis off;
end

%% get recon and ground truth
% load calibration
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
spectralPSF_3D = tifRead(filename);

% choose percentage of data to sample (or other sampling method here)
samplepercent = 0.1; %percent of pixels
[spectralPSF_2D,sampxy,sampfac] = randSamplePSF(spectralPSF_3D,samplepercent);
clear spectralPSF_3D

% load test spectrum
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0002.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';
spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);

% load wavelengths
load 'Datasets matFiles/calibrationFiles/calibrationWavelengths_fit.mat'

% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum_gt/partialSpectrum_2.spf2';
% fix the wavelength calibration on the ground truth spectrometer, resample
% to match calibrationWavelength_fit
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

sigma = 40
recon_gauss = gaussSVD(spectralPSF_2D,spectrumForRecon_sampled,sigma);

%% plot spectrum

figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend
box off
