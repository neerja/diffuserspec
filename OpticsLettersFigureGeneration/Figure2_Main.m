% Figure 2

% SETUP:
% Set home directory as BBOL server location:
%           Projects\smartOCT\Analysis\Data\CDMRPVRP\diffuserSpecData_OpticsLetters\organizedCodeData\Data
%
% The subFunctions folder within the above dir should be added to
% processing path.

% Since the code will ask you to select a main data
% folder and SSTM/SSTMbackground/calibrationFolder, it is not necessary to
% set the server folder as the main dir but it will reduce the number of clicks needed!

% Since data transfer from server is slow and data size of SSTM_3D is large
% (3.75GB) it is recommended to save the data to a local drive, but again
% this is not necessary, just for convenience 

%% Load data

disp('Select main data folder')
originalDir = uigetdir;

%load SSTM calibration matrix
disp('Select tape_SSTM_3D.tif')
    SSTM_fileName = uigetfile('*tif');
    SSTM = tifRead(SSTM_fileName);

%load background frame
disp('Select tape_SSTM_3D_background.tif')
    SSTM_bg_fileName = uigetfile('*tif');
    SSTM_bgFrame = tifRead(SSTM_bg_fileName);
    
[SSTM] = load_normalize_SSTM(SSTM,SSTM_bgFrame);

%% Load preCalibrated system data

% select calibrationFiles folder (located with
disp('Select calibrationFiles folder')
calibrationDir = uigetdir;
cd(calibrationDir)

load('calibrationOrig.mat');
load('calibrationCurve.mat');
load('signalRange.mat');
load('calibrationWavelengths_fit.mat');

cd(originalDir)

%%
% full mask sampling
    samplePercent = .5;
    innerRad = 0;
    outerRad = 0;

    [mask_full,maskCoordinates_full] = ...
        makeCircularMask(SSTM,innerRad,outerRad);
    
    [SSTM_full,samp_full] = applyMask_random...
        (SSTM,maskCoordinates_full,samplePercent);

   
% multipeak recon
    stepsize = ceil(344/6);
    offset = floor(round(344/6)-(round(344/6)/2));
    
    multiPeakVec = [stepsize-offset,(stepsize*2)-offset,(stepsize*3)-offset,...
        (stepsize*4)-offset,(stepsize*5)-offset,(stepsize*6)-offset];
    
    b_multiWave = SSTM_full(:,multiPeakVec);
    SSTM_full_downsamp = SSTM_full(:,1:2:end);
    
    sig = 1000;
        [recon_multipeak_full] = gaussSVD(SSTM_full_downsamp,...
        b_multiWave,sig);
    
%     recon_multipeak_full = recon_multipeak_full+(-min(recon_multipeak_full(:)));
    
%% Twopoint recon   

% shift spectral dimension so that middle wavelength is first
    SSTM_full_forSC = circshift(SSTM_full,172,2);

% calculate spectral correlation, limited to 50 wavelength shifts
    spectralCorrelation_full = ...
        calcSpectralCorrelation(SSTM_full_forSC, 50);
    
    SC_crop_full = mean(spectralCorrelation_full(:,1:20),1);
    SC_cropNorm_full = SC_crop_full/max(SC_crop_full);
    
    
    sig = 1000;
%downsample SSTM and select two wavelength frames for recon 
    
    b_full = mean(SSTM_full(:,[142 148]),2); 

    [recon_15pt_full] = gaussSVD(SSTM_full_downsamp,...
    b_full,sig);

%     recon_15pt_full = recon_15pt_full+(-min(recon_15pt_full));


%% plot

figure;
plot(SC_cropNorm_full)


figure;subplot(1,3,[1,2])
recon_multipeak_full_norm = recon_multipeak_full/max(recon_multipeak_full(:));
plot(calibrationWavelengths_fit(1:2:end),recon_multipeak_full_norm(1:end-1,:));
axis([-inf inf -.125 1.1])
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

subplot(1,3,3)
recon_15pt_full_norm = recon_15pt_full/max(recon_15pt_full(:));
plot(calibrationWavelengths_fit(132:2:160,:),recon_15pt_full_norm(66:80))
axis([-inf inf -.125 1.1])
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

