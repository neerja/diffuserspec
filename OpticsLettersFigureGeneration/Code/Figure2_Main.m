
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Figure2_Main.m
%
%  This script is meant to process data and generate figure content for 
%  Figure 2 of the DiffuserSpec manuscript
%
%  Other m-files required: 

%  Subfunctions: [tifRead.m, load_normalize_SSTM.m,
%  makeCircularMask.m, applyMask_random.m,gaussSVD.m]

%  MAT-files required: calibrationWavelengths_fit.mat

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

addpath('Code/subFunctions')

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

%%

    % The radius values are used to determine which region of the SSTM to
    % sample over. In this case, the 0's result in sampling the entire SSTM
    innerRad = 0;
    outerRad = 0;

    % makeCircularMask function uses the specified radii parameters to
    % create a binary mask used for SSTM sampling
    [mask_full,maskCoordinates_full] = ...
        makeCircularMask(SSTM,innerRad,outerRad);
    
    % Percentage of total SSTM to sample for the reconstruction, this value
    % is hard coded for 40000 pix.
    samplePercent = .18085*4;

    % applyMask_random function samples the SSTM according to the mask
    % generated above
    [SSTM_full,samp_full] = applyMask_random...
        (SSTM,maskCoordinates_full,samplePercent);

   
% Setup data for multipeak reconstruction

    % Automatically select 6 evenly spaced wavelength frames from the SSTM
    stepsize = ceil(344/6);
    offset = floor(round(344/6)-(round(344/6)/2));
    multiPeakVec = [stepsize-offset,(stepsize*2)-offset,(stepsize*3)-offset,...
        (stepsize*4)-offset,(stepsize*5)-offset,(stepsize*6)-offset];
    b_multiWave = SSTM_full(:,multiPeakVec);

    % downsample full SSTM from 344 frames to 172 frames, this is done so
    % that the multipeak data is not contained within the SSTM used for
    % reconstruction
    SSTM_full_downsamp = SSTM_full(:,1:2:end);
    
    % assign sigma filter value for gaussSVD reconstruction 
    sig = 250;

    % perform diffuse spectrum reconstruction of 6 peaks
    [recon_multipeak_full] = gaussSVD(SSTM_full_downsamp,...
        b_multiWave,sig);

% Twopoint recon       
    
    % select two wavelength frames spaced 2nm appart 
    b_full = mean(SSTM_full(:,[206 214]),2); 

    % perform diffuse spectrum reconstruction of two peaks
    [recon_15pt_full] = gaussSVD(SSTM_full_downsamp,...
    b_full,sig);

    % shift spectral dimension so that middle wavelength is first
    SSTM_full_forSC = circshift(SSTM_full,172,2);

%% Plot data

figure;subplot(1,3,[1,2])
recon_multipeak_full_norm = recon_multipeak_full./max(recon_multipeak_full);
recon_multipeak_full_norm(recon_multipeak_full_norm <0) = 0;
plot(calibrationWavelengths_fit(1:2:end),recon_multipeak_full_norm(1:end,:));hold on
stem(calibrationWavelengths_fit(multiPeakVec),ones(6,1),'k--')
axis([-inf inf 0 1.1])
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

subplot(1,3,3)
recon_15pt_full_norm = recon_15pt_full/max(recon_15pt_full(:));
recon_15pt_full_norm(recon_15pt_full_norm <0) = 0;
plot(calibrationWavelengths_fit(1:2:end,:),recon_15pt_full_norm);hold on
stem(calibrationWavelengths_fit([206 214]),ones(2,1),'k--')
axis([835 845 0 1.1])
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

%------------- END CODE --------------% 