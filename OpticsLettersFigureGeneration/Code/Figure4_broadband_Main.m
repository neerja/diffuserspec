
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Figure4_broadband_Main.m
%
%  This script is meant to process data and generate figure content for 
%  Figure 4 of the DiffuserSpec manuscript
%
%  Other m-files required: 

%  Subfunctions: [tifRead.m, load_normalize_SSTM.m,
%  makeCircularMask.m, applyMask_random.m, gaussSVD.m]

%  MAT-files required: calibrationWavelengths_fit.mat, calibrationCurve.mat
%  calibrationOrig.mat, signalRange.mat

%--------------------------------------------------------------------------
%
%  Author:          Joe Malone (joseph.d.malone@vanderbilt.edu)
%  Organization:    Bowden Biomedical Optics Lab, Vanderbilt University
%  Creation Date:   Sept. 2022
%  Last Modified:   2022/23/09
%
%  Matlab Style Template v 1.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%
% MAIN %
%%%%%%%%

%------------- BEGIN CODE --------------% 

% SETUP: PLEASE READ

% When code is run from the start, it will ask you to select a folder
% location. The folder to select is located on the Bowden lab server in the
% smartOCT projects folder (see directory location below). Once this
% directory is selected, the code will run to completion. Check command
% line during processing for status of code. 

% Projects\smartOCT\Publication\Articles\smartOCT_diffuseSpec\manuscript_Code_Data
 
%% Load data

clc

% Set working directory to server location that contains data and code for
% diffuserSpec manuscript
disp('Select folder specified in comments of this script (line 41)')
dataDir = uigetdir...
    ([],'Select folder specified in comments of this script (line 41)');
cd(dataDir)

addpath('Code\subFunctions')

cd('Data')

%load SSTM calibration matrix
disp('Loading SSTM')
    SSTM = tifRead('SSTM_3D.tif');

%load background frame
disp('Loading SSTM background frame')
    SSTM_bgFrame = tifRead('SSTM_background.tif');
    
% Normalize intensity of the SSTM calibration
disp('Normalizing Data - will take some time')
[SSTM] = load_normalize_SSTM(SSTM,SSTM_bgFrame);

%% Load preCalibrated system data

% Load spectrometer calibration data - this file provides information of
% the 'frame number to wavelength' calibration that was performed when the 
% system was first set up. This data is specific to the camera and
% custom monochrometer used in this system
load('calibrationFiles/calibrationWavelengths_fit.mat');

% ground truth spectrum collected on Thorlabs spectrometer 
load('broadbandData/groundTruth_wavelengths.mat');
load('broadbandData/groundTruth_spectrum.mat');


%% load broadband data

%load diffuse broadband image of Inphenix source
broadbandIm = tifRead('broadbandData/Inphenix_new_2.tif');

%load background frame
broadband_backgroundFrame = tifRead('broadbandData/superlum_background.tif');

%perform background subtraction
broadbandIm = abs(broadbandIm - broadband_backgroundFrame);

%% Broadband reconstruction 

    % Percentage of total SSTM to sample for the reconstruction, this value
    % is hard coded for 40000 pix.
    samplePercent = 5.35;

    % The radius values are used to determine which region of the SSTM to
    % sample over.
    innerRad = 1280;
    outerRad = 0;

    % makeCircularMask function uses the specified radii parameters to
    % create a binary mask used for SSTM sampling 
    [mask_corners,maskCoordinates_corners] = ...
        makeCircularMask(SSTM,innerRad,outerRad);

    % applyMask_random function samples the SSTM according to the mask
    % generated above
    [SSTM_corners,samp_corners] = applyMask_random...
        (SSTM,maskCoordinates_corners,samplePercent);

    % initialize variable for sampled broadbandIm matrix
    b_broadbandIm = zeros(size(SSTM_corners,1),1);

    % Sample broadband diffuse image using coordinates from above
    for coordinate = 1:size(SSTM_corners,1)
        b_broadbandIm(coordinate,:) = broadbandIm(maskCoordinates_corners...
            (samp_corners(coordinate,1),2),maskCoordinates_corners...
            (samp_corners(coordinate,1),1),:);
    end

    % assign sigma filter value and perform reconstruction
    sig =20;
    [recon_broadband] = gaussSVD(SSTM_corners(:,1:344),...
        b_broadbandIm,sig);

    % normalize to max value
    recon_broadband = recon_broadband/max(recon_broadband(10:330));

%% Plot data

% Fig. 4
figure;
plot(calibrationWavelengths_fit,recon_broadband,'b-.');hold on
plot(wavelength(41:553),spectrum(41:553),'k')
legend('DS','GT')
axis([786 870 -.1 1.2])
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')
