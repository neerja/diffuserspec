function [spectralPSF_3D] = import_spectralPSF_3D(filename,bgfilename,maxchannels)

% load calibration matrix
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D-001.tif'
spectralPSF_3D = tifRead(filename);

% load background frame
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/background_spectralPSF.tif'
backgroundFrame = tifRead(filename);

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

if narging>2
    spectralPSF_3D = spectralPSF_3D(:,:,1:maxchannels);
end

end