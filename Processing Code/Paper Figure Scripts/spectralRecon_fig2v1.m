%% spectralRecon_fig2v1.m
% Neerja Aggarwal & Joseph Malone
% July 13th, 2022
% Purpose: generate the recon for measurements.  

%% load calibration dataset (takes a minute)
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
spectralPSF_3D = tifRead(filename);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description

%% TODO:
% make a function that creates a sampling mask. takes in three parameters:
% example imagesize, min square,max square, number of sampling points. 
% use randSamplePSF then ignore any values inside min square. 

%% Create randomly sampled 2D spectralPSF
% choose percentage of data to sample
samplepercent = 0.1; %percent of pixels
[spectralPSF_2D,sampxy,sampfac] = randSamplePSF(spectralPSF_3D,samplepercent);

%% load midband spectrum
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0002.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';

spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);

% load wavelengths
load 'Datasets matFiles/calibrationFiles/calibrationWavelengths_fit.mat'

%% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum_gt/partialSpectrum_2.spf2';

% fix the wavelength calibration on the ground truth spectrometer, resample
% to match calibrationWavelength_fit
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

figure;
hold on
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt)
xlabel('Wavelength (nm)')
ylabel('normallized intensity')
legend('groundtruth')

%% set up gauss SVD parameters

sigma_range = [1:1:100]; % start and end values. 
errorvec = [];
tau = 1; % l1 error weight

% fig = figure;
for k1 = 1:length(sigma_range)
    sigma = sigma_range(k1)
    recon_gauss = gaussSVD(spectralPSF_2D,spectrumForRecon_sampled,sigma);
%     plot(calibrationWavelengths_fit,recon_gauss,'DisplayName',num2str(sigma))
%     hold on

% compute error
    l2error = norm(recon_gauss - partialSpectrum_gt);
    l1error = norm(recon_gauss,1);
    errorvec(k1) = l2error + tau*l1error;
end
% 
% plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
% xlabel('Wavelength (nm)'),
% ylabel('Intensity (arb unit)gma ')
% ylim([-0.2,1.2])
% xlim([782,868])
% legend

%% evaluate error and find min sigma value and plot
[val, sigind] = min(errorvec);
sigmaopt = sigma_range(sigind);


filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0002.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';
spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);
% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum_gt/partialSpectrum_2.spf2';
% fix the wavelength calibration on the ground truth spectrometer, resample
% to match calibrationWavelength_fit
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

recon_gauss = gaussSVD(spectralPSF_2D,spectrumForRecon_sampled,sigmaopt);
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend
%
figure; 
plot(sigma_range,errorvec);
xlabel('Sigma')
ylabel('Cost')
%% try the optimized value on partialspectrum #1

filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0001.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';

spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);

% load wavelengths
load 'Datasets matFiles/calibrationFiles/calibrationWavelengths_fit.mat'

% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum_gt/partialSpectrum_1.spf2';

% fix the wavelength calibration on the ground truth spectrometer, resample
% to match calibrationWavelength_fit
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

recon_gauss = gaussSVD(spectralPSF_2D,spectrumForRecon_sampled,sigmaopt);
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend

%% try the optimized value on partialspectrum #3

filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0003.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';

spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);

% load wavelengths
load 'Datasets matFiles/calibrationFiles/calibrationWavelengths_fit.mat'

% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum_gt/partialSpectrum_3.spf2';

% fix the wavelength calibration on the ground truth spectrometer, resample
% to match calibrationWavelength_fit
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

recon_gauss = gaussSVD(spectralPSF_2D,spectrumForRecon_sampled,sigmaopt);
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend

%% try on broadband spectrum 

filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadBand_Full/broadBand_#0003.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadBand_Full/broadband_background_#0004.tif';

spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);

% load wavelengths
load 'Datasets matFiles/calibrationFiles/calibrationWavelengths_fit.mat'

% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadBand_Full/broadBand3.spf2';

% fix the wavelength calibration on the ground truth spectrometer, resample
% to match calibrationWavelength_fit
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));


sigmaman = 20;
recon_gauss = gaussSVD(spectralPSF_2D,spectrumForRecon_sampled,sigmaman);
recon_gauss = recon_gauss./max(recon_gauss(2:end-2)); %ignore edge channels 
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend

%% try on broadband spectrum taken with broadest Pinhole
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadband_25nmBW/partialSpectrum_broadestPinhole_fullSpectrum.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadBand_Full/broadband_background_#0004.tif';

spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);

% load wavelengths
load 'Datasets matFiles/calibrationFiles/calibrationWavelengths_fit.mat'

% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadband_25nmBW/partialSpectrum_4_broadestPinhole.spf2';

% fix the wavelength calibration on the ground truth spectrometer, resample
% to match calibrationWavelength_fit
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt1 = partialSpectrum_gt/max(partialSpectrum_gt(:));

filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadband_25nmBW/partialSpectrum_4_broadestPinhole2.spf2';
% fix the wavelength calibration on the ground truth spectrometer, resample
% to match calibrationWavelength_fit
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(wavelengthOrig,wavelengthCorrected,signalRange,filename,calibrationWavelengths_fit);
partialSpectrum_gt2 = partialSpectrum_gt/max(partialSpectrum_gt(:));


sigmaman = 20;
recon_gauss = gaussSVD(spectralPSF_2D,spectrumForRecon_sampled,sigmaman);
recon_gauss = recon_gauss./max(recon_gauss(2:end-2)); %ignore edge channels 
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt1,'--','DisplayName','GT')
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt2,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend