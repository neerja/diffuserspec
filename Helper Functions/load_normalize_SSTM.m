% Purpose: subtract background and normalize the spectral transfer matrix 
% Inputs:
% spectralPSF_3D = 3D array containing calibration matrix with axes: x,y,lambda
% backgroundFrame = 2D image containg the background (dark) image
% Outputs:
% spectralPSF_3D = processed 3D array containing calibration matrix with axes: x,y,lambda
% normalizedIntensity = average intensity per wavelength (normalized)

function [spectralPSF_3D,normalizedIntensity] = load_normalize_SSTM(spectralPSF_3D,backgroundFrame)

%subtract background and resign negative values
    spectralPSF_3D = abs(spectralPSF_3D - backgroundFrame);

% get average intensity for each spectral PSF and normalize max to 1
    normalizedIntensity = mean(mean(spectralPSF_3D,1),2);
    normalizedIntensity = normalizedIntensity/max(normalizedIntensity);

% normalize spectral intensity and set max value to 1
    spectralPSF_3D = spectralPSF_3D./normalizedIntensity;
    spectralPSF_3D = spectralPSF_3D/max(spectralPSF_3D(:));

end
