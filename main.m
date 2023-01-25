% main.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: This script demonstrates basic principles presented in the
%  diffuserspec journal article
%
%  Instructions: You can run this entire script or section by section. Make
%  sure to set current folder to the diffuserspec/ directory. 
%
%  Subfunctions: All subfunctions are located in ./HelperFunctions/
%  applyMask_random.m, calcSpectralCorrelation.m, gaussSVD.m,
%  load_normalize_SSTM.m, makeCircularMask.m, tifRead.m
%
%  Data files required: All data files are located in ./Data/ 
%  Code and data can be downloaded from https://osf.io/b79an/
%--------------------------------------------------------------------------
%
%  Author:          Joe Malone (joseph.d.malone@vanderbilt.edu, Neerja
%                   Aggarwal (neerja@berkeley.edu)
%  Organization:    Bowden Biomedical Optics Lab (Vanderbilt University) &
%                   Computational Imaging Lab (University of California -
%                   Berkeley)
%  Creation Date:   2022/20/12
%  Last Modified:   2022/20/12
%  Cite as:         J. Malone, N. Aggarwal, L. Waller, A. Bowden.
%                   "DiffuserSpec: spectroscopy with Scotch tape," 
%                   Optics Letters (2022)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%
% MAIN %
%%%%%%%%

%------------- BEGIN CODE --------------% 

%% 1. Load calibration files

% The calibration data is located in the Data folder.  

% load system calibration files
load('./Data/calibrationFiles/calibrationWavelengths_fit.mat');

% assign file locations
SSTM_fileName = './Data/SSTM_3D.tif';
SSTM_background_fileName = './Data/SSTM_background.tif';

%load SSTM calibration matrix
disp('Loading SSTM...')
    SSTM = tifRead(SSTM_fileName);

%load background frame
disp('Loading SSTM background frame...')
    SSTM_backgroundFrame = tifRead(SSTM_background_fileName);
    
% Normalize intensity of the SSTM calibration
disp('Normalizing Data - will take some time')
    [SSTM] = load_normalize_SSTM(SSTM,SSTM_backgroundFrame);

%% 2. Load measurement

% ground truth spectrum collected on Thorlabs spectrometer 
load('./Data/broadbandMeasurements/Inphenix_groundTruth_wavelengths.mat');
load('./Data/broadbandMeasurements/Inphenix_groundTruth_spectrum.mat');

% crop relevant wavelength range (786-870nm)
groundTruth_wavelengths = groundTruth_wavelengths(41:553);
groundTruth_spectrum = groundTruth_spectrum(41:553);

% assign filenames for measurement data
measurement_fileName = './Data/broadbandMeasurements/Inphenix_1.tif';
measurementBackground_fileName = './Data/broadbandMeasurements/Inphenix_background.tif';

% load diffuse broadband image of Inphenix source
broadbandMeasurement = tifRead(measurement_fileName);

% load background frame
broadbandMeasurement_background = tifRead(measurementBackground_fileName);

% perform background subtraction
broadbandMeasurement = abs(broadbandMeasurement - broadbandMeasurement_background);

% visualize speckle measurement, with intensity range 0:3E4
figure('Renderer', 'painters', 'Position', [100 100 1600 600])
subplot(1,2,1)
imagesc(SSTM(:,:,129));colormap turbo; colorbar
title('DiffuserSpec monochromatic (820nm) speckle measurment')
xlabel('X-pixel')
ylabel('Y-pixel')

subplot(1,2,2)
imagesc(broadbandMeasurement,[500 3E4]);colormap turbo; colorbar
title('DiffuserSpec broadband speckle measurment')
xlabel('X-pixel')
ylabel('Y-pixel')


%% 3. Subsampling

% Percentage of total SSTM to sample for the reconstruction, this value
% is hard coded for 40000 pix.
samplePercent = 5.35;

% The radius values are used to determine which region of the SSTM to
% sample over.
innerRad = 1280;
outerRad = 0;

% makeCircularMask function uses the specified radii parameters to
% create a binary mask used for SSTM sampling 
[mask_corners,maskCoordinates_outerRad] = ...
    makeCircularMask(SSTM,innerRad,outerRad);

% applyMask_random function samples the SSTM according to the mask
% generated above
[SSTM_outerRad,samp_corners] = applyMask_random...
    (SSTM,maskCoordinates_outerRad,samplePercent);

% initialize variable for sampled broadbandIm matrix
broadbandMeasurement_sampled = zeros(size(SSTM_outerRad,1),1);

% Sample broadband diffuse image using coordinates from above
for coordinate = 1:size(SSTM_outerRad,1)
    broadbandMeasurement_sampled(coordinate,:) = broadbandMeasurement(maskCoordinates_outerRad...
        (samp_corners(coordinate,1),2),maskCoordinates_outerRad...
        (samp_corners(coordinate,1),1),:);
end


%% 4. Run algorithm to obtain reconstructed spectrum

% assign sigma filter value and perform reconstruction, suggested 30
sig =20;

% use the gaussSVD algorithm (pseudoinverse) to estimate the spectrum
[recon_broadband] = gaussSVD(SSTM_outerRad,...
    broadbandMeasurement_sampled,sig);

% normalize to max value
recon_broadband = recon_broadband/max(recon_broadband);

% plot reconstructed and ground truth spectrum
figure;
plot(calibrationWavelengths_fit,recon_broadband,'b-.');hold on
plot(groundTruth_wavelengths,groundTruth_spectrum./max(groundTruth_spectrum),'k')
legend('DiffuserSpec','Ground Truth')
axis([786 870 -.1 1.2])
title('DiffuserSpec reconstruction of Inphenix SLED source')
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

%% 5. analyze the spectral correlation for this diffuser matrix

% assign spectral shift variable based on calibrated wavelength values and
% use index 172 as central point to mimic spectral correlation calculation
% above
spectralShift = calibrationWavelengths_fit-calibrationWavelengths_fit(172);
spectralShift = circshift(spectralShift,173);

% shift spectral dimension so that middle wavelength is first
SSTM_outerRad_forSC = circshift(SSTM_outerRad,172,2);

% calculate spectral correlation, 
spectralCorrelation_outerRad = ...
    calcSpectralCorrelation(SSTM_outerRad_forSC, 1,344);

SC_crop_outerRad = mean(spectralCorrelation_outerRad(:,1:25),1);
SC_cropNorm_outerRad = SC_crop_outerRad/max(SC_crop_outerRad);

figure('Renderer', 'painters', 'Position', [100 700 1600 600])
subplot(1,2,1)
imagesc(calibrationWavelengths_fit,1:size(SSTM_outerRad,1),SSTM_outerRad);
colormap turbo; colorbar
title('Sampled SSTM')
ylabel('Pixel Number')
xlabel('Wavelength (nm)')

subplot(1,2,2)
plot(spectralShift(1:25),SC_cropNorm_outerRad);
axis([0 5 -.1 1])
title('Spectral correlation (C) of sampled SSTM')
xlabel('Spectral Shift (nm)')
ylabel('Norm. Spectral Correlation (C)')
    
    
%------------- END CODE --------------% 
