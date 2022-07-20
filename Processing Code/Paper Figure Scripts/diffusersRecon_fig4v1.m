%% diffuserRecon_fig4v1.m
% Neerja Aggarwal & Joseph Malone
% July 13th, 2022
% Purpose: generate the recon for measurements.  

%% load calibration dataset for 120grit (takes a minute)
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
spectralPSF_3D = tifRead(filename);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description

samplepercent = 0.1;
maskbox = 500;

% randomly sample the PSF with a masked box in center 
[spectralPSFrand, samp_ind_nonan] = randSamplePSF_mask(spectralPSF_3D,samplepercent,maskbox);

%% load midband spectrum
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadBand_Full/broadband_#0001.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadBand_Full/broadband_background_#0004.tif';

spectrumForRecon = loadTestMeasurement(filename,bgfilename);
spectrumForRecon_sampled = spectrumForRecon(samp_ind_nonan);

% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/broadBand_Full/broadBand1.spf2';

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

% recon the spectrum
sigma = 20;
recon_gauss = gaussSVD(spectralPSFrand,spectrumForRecon_sampled,sigma);
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend

%% try on TAPE
% load calibration dataset for tape (takes a minute)
filename = './Raw Data/Tape/tape_spectralTM_3D-001.tif'
bgfilename = './Raw Data/Tape/tape_spectralTM_background.tif'

spectralPSF_3D = import_spectralPSF_3D(filename,bgfilename);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description

samplepercent = 0.1;
maskbox = 500;

% randomly sample the PSF with a masked box in center 
[spectralPSFrand, samp_ind_nonan] = randSamplePSF_mask(spectralPSF_3D,samplepercent,maskbox);

%% load and recon broadband spectrum
filename = './Raw Data/Tape/broadband/broadband_full_#0001.tif';
bgfilename = './Raw Data/Tape/broadband/broadband_full_background_#0006.tif';

spectrumForRecon = loadTestMeasurement(filename,bgfilename);
spectrumForRecon_sampled = spectrumForRecon(samp_ind_nonan);

% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/Tape/broadband_gt/broadband_source2_gt.spf2';

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

% recon the spectrum
sigma = 10;
recon_gauss = gaussSVD(spectralPSFrand,spectrumForRecon_sampled,sigma);
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend

%% try with DOE

% load calibration dataset for tape (takes a minute)
filename = './Raw Data/DOE/DOE_spectralTM_3D-002.tif'
bgfilename = './Raw Data/DOE/DOE_spectralTM_background.tif'

spectralPSF_3D = import_spectralPSF_3D(filename,bgfilename);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description

samplepercent = 0.1;
maskbox = 500;

% randomly sample the PSF with a masked box in center 
[spectralPSFrand, samp_ind_nonan] = randSamplePSF_mask(spectralPSF_3D,samplepercent,maskbox);

%% load and recon broadband spectrum
filename = './Raw Data/DOE/broadBand/broadband_source1_#0003.tif';
bgfilename = './Raw Data/DOE/broadBand/broadband_background_source1_#0006.tif';

spectrumForRecon = loadTestMeasurement(filename,bgfilename);
spectrumForRecon_sampled = spectrumForRecon(samp_ind_nonan);

% load gt spectrum 
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/DOE/broadband_gt/broadband_source1_gt.spf2';

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

%% recon the spectrum
sigma = 20;
recon_gauss = gaussSVD(spectralPSFrand,spectrumForRecon_sampled,sigma);
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend

%% try find a single measurement.

sigma = 40;
recon_gauss = gaussSVD(spectralPSFrand,spectralPSFrand(:,200),sigma);
figure;
plot(calibrationWavelengths_fit,recon_gauss,'DisplayName','Recon'); hold on;
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'--','DisplayName','GT')
xlabel('Wavelength (nm)'),
ylabel('Intensity (arb unit)')
ylim([-0.2,1.2])
xlim([782,868])
legend


