%% load diffuser spec data

%load calibration matrix
disp('Select calibration matrix file')
spectrumForRecon = uigetfile('.tif');
spectralPSF_3D = tifRead(spectrumForRecon);

%% load background frame
disp('Select background file')
backgroundFilename = uigetfile('.tif');
backgroundFrame = tifRead(backgroundFilename);

%subtract background and resign negative values
spectralPSF_3D = abs(spectralPSF_3D - backgroundFrame);

%% Spectral normalization of calibration data

% get average intensity for each spectral PSF and set range between [0,1]
normalizedIntensity = mean(mean(spectralPSF_3D,1),2);
normalizedIntensity = normalizedIntensity/max(normalizedIntensity);

% normalize spectral intensity and set max value to 1
spectralPSF_3D = spectralPSF_3D./normalizedIntensity;
spectralPSF_3D = spectralPSF_3D/max(spectralPSF_3D(:));

%% process ground truth wavelength data 

load('calibrationOrig.mat');
load('calibrationCurve.mat');
load('signalRange.mat');

%read in wavelength data
[spectrumData_gt,wavelengthData_gt] = readSPF2_sequence(calibrationCurve,calibrationOrig,signalRange);

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

%% alternative wavelength data (Neerja added)
load ./JoeData_07062022/calibrationWavelengths.mat
load ./JoeData_07062022/calibrationWavelengths_polyfit.mat
calibrationWavelengths_fit = calibrationWavelengths_polyfit;

%% Create randomly sampled 2D spectralPSF

% choose percentage of data to sample
samplepercent = .1;

[N1,N2,N3] = size(spectralPSF_3D);
sampfac = round(sqrt(N1*N2*(samplepercent/100)));

samp_xy = zeros(sampfac,2);
samp_xy(:,1) = ceil(N1*rand(sampfac,1));
samp_xy(:,2) = ceil(N2*rand(sampfac,1));

samp_full = spectralPSF_3D(samp_xy(:,1),samp_xy(:,2),:);
samp_line = reshape(samp_full,sampfac^2,N3)';

%% spectral correlation calc

% spectral correlation calculation (taken from Wang, Menon OE 2014)
    midLine = round(size(samp_line,1)/2);
    PSFt = mean(samp_line(midLine,:).*samp_line,2);
    PSFb_resize = mean(samp_line(midLine,:),2).*mean(samp_line,2);

    spectralCorrelation_resize = (PSFt./PSFb_resize)-1;

    spectralCorrelation_resizeNormalized = ...
        spectralCorrelation_resize/max(spectralCorrelation_resize);
    
    figure;plot(calibrationWavelengths_fit(:,173:end)',...
        spectralCorrelation_resizeNormalized(173:344,:));
    axis([-inf inf 0 1])
    xlabel('Wavelength (nm)')
    ylabel('Spectral Correlation')

%% Reconstruct partial spectrum

disp('Select spectrum for recon file')
spectrumForReconFile = uigetfile('.tif');
spectrumForRecon = tifRead(spectrumForReconFile);

disp('Select background spectrum for recon file')
bgSpectrumForReconFile = uigetfile('.tif');
bgSpectrumForRecon = tifRead(bgSpectrumForReconFile);

spectrumForRecon = abs(spectrumForRecon - bgSpectrumForRecon);
spectrumForRecon_sampled = spectrumForRecon(samp_xy(:,1),samp_xy(:,2));
spectrumForRecon_sampled = reshape(spectrumForRecon_sampled,sampfac^2,1);

thresh = 0.2; %2 is optimal
[recon,Alr,Ainv,Dprime,Dprime_thresh] = TruncSVD(samp_line',...
    spectrumForRecon_sampled,thresh);
recon = abs(recon/max(recon(:)));

%load groundtruth spectrum
% [partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
%     readSPF2_withCalibration(calibrationCurve,calibrationOrig,signalRange);

% partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

% partialSpectrum_gt = imresize(partialSpectrum_gt,...
%     [344 size(partialSpectrum_gt,2)]);
% partialSpectrumWavelengths_gt = imresize(partialSpectrumWavelengths_gt,...
%     [344 size(partialSpectrumWavelengths_gt,2)]);

figure;plot(calibrationWavelengths_fit',recon(1:344,:));hold on
xlabel('Wavelength (nm)')
ylabel('normallized intensity')
legend('Recon','groundtruth')



%% NOT WORKING: (missing files)
% load groundtruth spectrum
disp('Select ground truth for measurement')
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(calibrationCurve,calibrationOrig,signalRange);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

partialSpectrum_gt = imresize(partialSpectrum_gt,...
    [344 size(partialSpectrum_gt,2)]);
partialSpectrumWavelengths_gt = imresize(partialSpectrumWavelengths_gt,...
    [344 size(partialSpectrumWavelengths_gt,2)]);

plot(partialSpectrumWavelengths_gt,partialSpectrum_gt)
xlabel('Wavelength (nm)')
ylabel('normallized intensity')
legend('groundtruth')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% sub-functions