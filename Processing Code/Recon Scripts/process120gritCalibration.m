%% process120gritCalibration.m
% Neerja Aggarwal
% June 28th, 2022
% Purpose: make spectralPSF_3D_norm.mat file containing normalized
% calibration matrix and wavelengths. 

% IMPORTANT - run this script with the repo base directory as the current 
% folder: ex: diffuserspec/ 

%% load calibration data (takes a minute)
% load calibration matrix
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D-001.tif'
spectralPSF_3D = tifRead(filename);

% load background frame
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/background_spectralPSF.tif'
backgroundFrame = tifRead(filename);
%% process calibration dataset

%subtract background and resign negative values
spectralPSF_3D = abs(spectralPSF_3D - backgroundFrame);
% get average intensity for each spectral PSF and set range between [0,1]
normalizedIntensity = mean(mean(spectralPSF_3D,1),2);
normalizedIntensity = normalizedIntensity/max(normalizedIntensity);

% normalize spectral intensity so they are roughly the same per wavelength
spectralPSF_3D = spectralPSF_3D./normalizedIntensity;
% normalize calibration data set so max pixel value is one. (To avoid low
% intensity reconstructed spectra)
spectralPSF_3D = spectralPSF_3D/max(spectralPSF_3D(:));

%% load calibration wavelengths
load('./Datasets matFiles/120grit/calibrationOrig.mat');
load('./Datasets matFiles/120grit/calibrationCurve.mat');
load('./Datasets matFiles/120grit/signalRange.mat');

%read in wavelength data
directory = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationMatrix_gt/'
[spectrumData_gt,wavelengthData_gt] = readSPF2_sequence(calibrationCurve,calibrationOrig,signalRange,directory);

%find location of wavelength peaks
[~,I] = max(spectrumData_gt,[],1);

%find wavelengths corresponding to peak index
calibrationWavelengths = wavelengthData_gt(I);

% fit 2nd order polynomial to wavelength peaks to smooth data and
% compensate for max peak identification errors
x = 1:size(calibrationWavelengths,1);
fit = polyfit(1:size(calibrationWavelengths,1),...
    calibrationWavelengths,2);
calibrationWavelengths_fit = polyval(fit,x);

%% save processed spectralPSF_120grit.mat
% remove extra calibration wavelengths and match the gt number:
[~,N] = size(calibrationWavelengths_fit);
spectralPSF_3D = spectralPSF_3D(:,:,1:N);
saveFilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif';
description = 'normalized unsampled calibration dataset with fitted wavelengths for 120grit_2022-06-16'
for k1 = 1:N
    img = squeeze(spectralPSF_3D(:,:,k1));
    k1
    if k1 == 1
        % First slice:
        imwrite(img,saveFilename,'TIFF');
    else
        imwrite(img,saveFilename,'TIFF','WriteMode','append');
    end
end

saveFilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat';
save(saveFilename, 'calibrationWavelengths_fit','description')