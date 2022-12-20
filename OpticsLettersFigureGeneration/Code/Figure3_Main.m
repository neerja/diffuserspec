
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Figure3_Main.m
%
%  This script is meant to process data and generate figure content for 
%  Figure 3 of the DiffuserSpec manuscript
%
%  Other m-files required: 

%  Subfunctions: [tifRead.m, load_normalize_SSTM.m,
%  makeCircularMask.m, applyMask_random.m,gaussSVD.m, calcSpectralCorrelation.m]

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

%% Figure 4a: calculate 2D spatial correlation data

% initialize spectral correlation variable
    spectralCorrelation = zeros(2160,2560,50);

for row = 1:2560

    % shift spectral dimension so that middle wavelength is first
    SSTM_insideMiddleCircle = circshift(squeeze(SSTM(:,row,:)),172,2);

    for shift = 1:50

        % calculate spectral correlation numerator and denomenator (see
        % Menon paper for full expression)
        psftop = mean((SSTM_insideMiddleCircle.*...
            circshift(SSTM_insideMiddleCircle,shift,2)),2);
        psfbottom = mean(SSTM_insideMiddleCircle,2).*...
            mean(circshift(SSTM_insideMiddleCircle,shift,2),2);
    
        spectralCorrelation(:,row,shift)= (psftop./psfbottom) - 1;
    end
    
    % Display status of calculations 
    if mod(row,100) == 0
        disp([num2str(row),'/2560'])
    else
    end

end

% normalize correlation data to max of 1
    spectralCorrelation_norm = spectralCorrelation./max(spectralCorrelation,[],3);

%set threshold and wavelength maping 
    thresh = 0.5;
    
% assign spectral shift variable based on calibrated wavelength values and
% use index 172 as central point to mimic spectral correlation calculation
% above
    spectralShift = calibrationWavelengths_fit-calibrationWavelengths_fit(172);
    spectralShift = circshift(spectralShift,173);

%initialize spectral resolution variable
    spectralRes = zeros(2160,2560);

for colNum = 1:size(spectralCorrelation_norm,1)
    for rowNum = 1:size(spectralCorrelation_norm,2)

        % find index of values for each pixel of 2D spectral correlation
        % that is less than the threshold 
        [~,ind] = find(spectralCorrelation_norm(colNum,rowNum,:)<thresh);
    
        % find resolution shift that corresponds to the first index value
        if isempty(ind)
            %if no index, resolution set to max shift
            spectralRes(colNum,rowNum) = spectralShift(50);
        else
            spectralRes(colNum,rowNum) = spectralShift(ind(1));
        end
    end
    
    % Display status of calculations
    if mod(colNum,100) == 0
        disp([num2str(colNum),'/2160'])
    else
    end
    
end

% apply gaussian filter to smooth data
    spectralRes_gaussFilt = imgaussfilt(spectralRes,50);

%% Mask sampling spectral correlation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Inside of middle circle region %%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Percentage of total SSTM to sample for the reconstruction, this value
    % is hard coded for 40000 pix.
    samplePercent = 1.27335*4;

    % The radius values are used to determine which region of the SSTM to
    % sample over.
    innerRad = 0;
    outerRad = 500;

    % makeCircularMask function uses the specified radii parameters to
    % create a binary mask used for SSTM sampling
    [mask_insideMiddleCircle,maskCoordinates_insideMiddleCircle] = ...
        makeCircularMask(SSTM,innerRad,outerRad);

    % applyMask_random function samples the SSTM according to the mask
    % generated above
    [SSTM_insideMiddleCircle,samp_insideMiddleCircle] = applyMask_random...
        (SSTM,maskCoordinates_insideMiddleCircle,samplePercent);
    
    % shift spectral dimension so that middle wavelength is first
    SSTM_insideMiddleCircle_forSC = circshift(SSTM_insideMiddleCircle,172,2);

    % calculate spectral correlation, limited to 50 wavelength shifts
    spectralCorrelation_insideMiddleCircle = ...
        calcSpectralCorrelation(SSTM_insideMiddleCircle_forSC, 50);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Outside of middle Box region %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    samplePercent = .2108*4;
    innerRad = 500;
    outerRad = 0;
    [mask_outsideMiddleCircle,maskCoordinates_outsideMiddleCircle] = ...
        makeCircularMask(SSTM,innerRad,outerRad);
    [SSTM_outsideMiddleCircle,samp_outsideMiddleCircle] = applyMask_random...
        (SSTM,maskCoordinates_outsideMiddleCircle,samplePercent);
    
    % shift spectral dimension so that middle wavelength is first
    SSTM_outsideMiddleCircle_forSC = circshift(SSTM_outsideMiddleCircle,172,2);

    % calculate spectral correlation, limited to 50 wavelength shifts
    spectralCorrelation_outsideMiddleCircle = ...
        calcSpectralCorrelation(SSTM_outsideMiddleCircle_forSC, 50);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          Corners             %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    samplePercent = 1.3245*4;
    innerRad = 1280;
    outerRad = 0;
    [mask_corners,maskCoordinates_corners] = ...
        makeCircularMask(SSTM,innerRad,outerRad);
    [SSTM_corners,samp_corners] = applyMask_random...
        (SSTM,maskCoordinates_corners,samplePercent);
    
    % shift spectral dimension so that middle wavelength is first
    SSTM_corners_forSC = circshift(SSTM_corners,172,2);

    % calculate spectral correlation, limited to 50 wavelength shifts
    spectralCorrelation_corners = ...
        calcSpectralCorrelation(SSTM_corners_forSC, 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Average along pixel locations and normalize each correlation function to 1
    SC_crop_insideMiddleCircle = mean(spectralCorrelation_insideMiddleCircle(:,1:20),1);
    SC_cropNorm_insideMiddleCircle = SC_crop_insideMiddleCircle/max(SC_crop_insideMiddleCircle);

    SC_crop_outsideMiddleCircle = mean(spectralCorrelation_outsideMiddleCircle(:,1:20),1);
    SC_cropNorm_outsideMiddleCircle = SC_crop_outsideMiddleCircle/max(SC_crop_outsideMiddleCircle);

    SC_crop_corners = mean(spectralCorrelation_corners(:,1:20),1);
    SC_cropNorm_corners = SC_crop_corners/max(SC_crop_corners);

    
%%  point recons

%select sigma value for recon algoithm
    sig =250;
    
% inside middleBox
    %downsample SSTM and select two wavelength frames for recon 
    SSTM_insideMiddleCircle_downsamp = SSTM_insideMiddleCircle(:,1:2:end);
    b_insideMiddleCircle = mean(SSTM_insideMiddleCircle(:,[142]),2); 

    % perform reconstruction
    [recon_15pt_insideMiddleCircle] = gaussSVD(SSTM_insideMiddleCircle_downsamp,...
    b_insideMiddleCircle,sig);
    % soft threshold neg values
    recon_15pt_insideMiddleCircle(recon_15pt_insideMiddleCircle<0) = 0;
    
% outsidemiddleBox
    SSTM_outsideMiddleCircle_downsamp = SSTM_outsideMiddleCircle(:,1:2:end);
    res15nm_outsideMiddleCircle = mean(SSTM_outsideMiddleCircle(:,[142]),2); 
    
    [recon_15pt_outsideMiddleCircle] = gaussSVD(SSTM_outsideMiddleCircle_downsamp,...
    res15nm_outsideMiddleCircle,sig);
    recon_15pt_outsideMiddleCircle(recon_15pt_outsideMiddleCircle<0) = 0;

% corner
    SSTM_corners_downsamp = SSTM_corners(:,1:2:end);
    res15nm_corners = mean(SSTM_corners(:,[142]),2); 
    
    [recon_15pt_corners] = gaussSVD(SSTM_corners_downsamp,...
    res15nm_corners,sig);
    recon_15pt_corners(recon_15pt_corners<0) = 0;

    % Normalize reconstructions to max value of 1
    recon_15pt_insideMiddleCircle = recon_15pt_insideMiddleCircle/max(recon_15pt_insideMiddleCircle(:));
    recon_15pt_outsideMiddleCircle = recon_15pt_outsideMiddleCircle/max(recon_15pt_outsideMiddleCircle(:));
    recon_15pt_corners = recon_15pt_corners/max(recon_15pt_corners(:));
    
    
%% Fig. 3 plots

% Fig. 3(a)
figure;
subplot(2,2,1);
imagesc(spectralRes_gaussFilt); colormap turbo
    
% Fig. 3(c)
subplot(2,2,3)
    plot(spectralShift(1:20),SC_cropNorm_insideMiddleCircle);hold on
    plot(spectralShift(1:20),SC_cropNorm_outsideMiddleCircle);
    plot(spectralShift(1:20),SC_cropNorm_corners);
    plot(spectralShift(1:20),ones(20,1)/2)
    axis([0 5 -.1 1])
    xlabel('Spectral Shift (nm)')
    ylabel('Norm. Spectral Correlation (C)')
    legend('Inner Circle','Outer Circle','Corners')

% Fig. 3(d)
subplot(2,2,4)
    figure;
    plot(calibrationWavelengths_fit(1:2:end,:),...
        recon_15pt_insideMiddleCircle);hold on
    plot(calibrationWavelengths_fit(1:2:end,:),...
        recon_15pt_outsideMiddleCircle);
    plot(calibrationWavelengths_fit(1:2:end,:),...
        recon_15pt_corners);
    axis([815.5 831.5 0 1.1])
    legend('Inner Circle','Outer Circle','Corners')
    xlabel('Wavelength (nm)')
    ylabel('Norm. Intensity')
    
% visualize masks for Fig. 3(b)
figure;
    subplot(2,2,1)
    imagesc(mask_insideMiddleCircle);colormap gray
    subplot(2,2,2)
    imagesc(mask_outsideMiddleCircle);colormap gray
    subplot(2,2,3)
    imagesc(mask_corners);colormap gray

    %------------- END CODE --------------% 