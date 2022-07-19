%% radialProcess.m
% Neerja Aggarwal
% July 18th, 2022

% Use the radial transformed data to do broadband recon. 
% 1.  Import the 3D matrix
% 2.  Import the measurement
% 3.  Import the ground truth
% 4.  Use gauss SVD and try the recon with sigma = 20 or 40

numWavelengths = 344; 

%% load calibration dataset (takes a minute)
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/radialData/radialReslice_120grit.tif'
spectralPSF_3D = tifRead(filename);

%subtract background and resign negative values
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/radialData/radialReslice_120grit_background.tif'
spectralPSF_bg = tifRead(bgfilename);
spectralPSF_3D = abs(spectralPSF_3D-spectralPSF_bg);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description

% remove the last frame cause it's empty
spectralPSF_3D = spectralPSF_3D(:,:,1:numWavelengths);
[N1,N2,N3] = size(spectralPSF_3D);

% get average intensity for each spectral PSF and set range between [0,1]
normalizedIntensity = mean(mean(spectralPSF_3D,1),2);
normalizedIntensity = normalizedIntensity/max(normalizedIntensity);

% normalize spectral intensity so they are roughly the same per wavelength
spectralPSF_3D = spectralPSF_3D./normalizedIntensity;
% normalize calibration data set so max pixel value is one. (To avoid low
% intensity reconstructed spectra)
spectralPSF_3D = spectralPSF_3D/max(spectralPSF_3D(:));

clear spectralPSF_bg

%% load test spectrum 

% % load partial broadband
% filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/radialData/RadialReslice_120grit_partialBroadband.tif';
% bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/radialData/RadialReslice_120grit_partialBroadband_background.tif';
% measurement = tifRead(filename);
% bg = tifRead(bgfilename);
% measurement = measurement-bg;
% clear bg

% load full broadband
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/radialData/RadialReslice_120grit_broadband1.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/radialData/RadialReslice_120grit_broadband1_background.tif';
measurement = mean(tifRead(filename),3);
bg = tifRead(bgfilename);
bg(:,N2) = bg(:,end); % copy over last row of bg to make it the right size.
measurement = measurement-bg;
clear bg

%% load gt for full broadband
load './Datasets matFiles/calibrationFiles/wavelength_gt.mat'
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/radialData/broadBand1.spf2';

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

%% visualize the calibration matrices

N1max = 1000;
% for k1 = 1:N1max

for k1 = [10,100,500,1000]
    A = squeeze(spectralPSF_3D(:,k1,:)); % only 360 samples. 
    figure
    imagesc(A)
    colorbar()
    xlabel('Wavelength Channel')
    ylabel('Pixel')
    name = horzcat('Radius = ', num2str(k1));
    title(name)
end

%% try  gauss SVD recon 
rmin = 100;
rmax = 1000;
rstep = 100;
rvec = [rmin:rstep:rmax];
sigma = 2;
figure;
hold on;
for k = 1:length(rvec)
    r = rvec(k)
    A = squeeze(spectralPSF_3D(:,r,:));
    b = measurement(:,r);
    x = gaussSVD(A,b,sigma);
    plot(calibrationWavelengths_fit,x)
end

%% try ridge regression
rmin = 100;
rmax = 1000;
rstep = 100;
rvec = [rmin:rstep:rmax];
lambda = 1;
figure;
hold on;
for k = 1:length(rvec)
    r = rvec(k)
    A = squeeze(spectralPSF_3D(:,r,:));
    b = measurement(:,r);
    x = gaussSVD(A,b,lambda);
    plot(calibrationWavelengths_fit,x)
end

figure;
hold on
plot(partialSpectrumWavelengths_gt,partialSpectrum_gt)
xlabel('Wavelength (nm)')
ylabel('normallized intensity')
legend('groundtruth')

%% Conclusion:

% Radial sampling greatly reduces the number of measurements (ex: 360 theta
% values).  Might be too few to be able to properly do recon.  Can try
% doing a band? 

%%  try banded radial (kind of works at specific radii)
rmin = 500;
rmax = 1000;
rstep = 100;
rvec = [rmin:rstep:rmax];
rband = 10; %how many radii to use
sigma = 20;
figure;
hold on;

[M1,M2] = size(measurement);

for k = 1:length(rvec)
    %choose radial band
    r = rvec(k)
    A = spectralPSF_3D(:,r:r+rband-1,:);
    b = measurement(:,r:r+rband-1);
    % turn calibration matrix into 2D
    Ap = reshape(A,[M1*rband,N3]);
    % vectorize the measurement
    bp = b(:);
    % do recon
    x = gaussSVD(Ap,bp,sigma);
    plot(calibrationWavelengths_fit,x,'DisplayName',num2str(r))
    ylim([0,1]);
end

plot(partialSpectrumWavelengths_gt,partialSpectrum_gt,'k--','DisplayName','ground truth')
xlabel('Wavelength (nm)')
ylabel('normallized intensity')
legend()
title('120 grit broadband at different radius')