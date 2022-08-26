% DiffuserSpec Optics Letters manuscript processing

% FIGURE 4

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

%% Figure 4a: calculate 2D spatial correlation data

% initialize variables
    spectralCorrelation = zeros(2160,2560,50);

for row = 1:2560

    % shift spectral dimension so that middle wavelength is first
    SSTM_insideMiddleCircle = circshift(squeeze(SSTM(:,row,:)),172,2);

    for shift = 1:50

        % calculate spectral correlation numerator and denomenator (see
        % Menon paper for full expression)
        psftop = mean((SSTM_insideMiddleCircle.*circshift(SSTM_insideMiddleCircle,shift,2)),2);
        psfbottom = mean(SSTM_insideMiddleCircle,2).*mean(circshift(SSTM_insideMiddleCircle,shift,2),2);
    
        spectralCorrelation(:,row,shift)= (psftop./psfbottom) - 1;
    end
    
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
    
    % display for timing
    if mod(colNum,100) == 0
        disp([num2str(colNum),'/2160'])
    else
    end
    
end

%apply gaussian filter to smooth data
    spectralRes_gaussFilt = imgaussfilt(spectralRes,50);

%% Mask sampling spectral correlation

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Inside of middle circle region %%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    samplePercent = 1.28;
    innerRad = 0;
    outerRad = 500;
    [mask_insideMiddleCircle,maskCoordinates_insideMiddleCircle] = ...
        makeCircularMask(SSTM,innerRad,outerRad);
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

    samplePercent = .215;
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
%%%%%%            Corners           %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    samplePercent = 1.33;
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

    
%%  Twopoint recons

%select sigma value for recon algoithm
    sig = 1000;
    
% inside middleBox
    %downsample SSTM and select two wavelength frames for recon 
    SSTM_insideMiddleCircle_downsamp = SSTM_insideMiddleCircle(:,1:2:end);
    b_insideMiddleCircle = mean(SSTM_insideMiddleCircle(:,[142 148]),2); 

    [recon_15pt_insideMiddleCircle] = gaussSVD(SSTM_insideMiddleCircle_downsamp,...
    b_insideMiddleCircle,sig);
    
% outsidemiddleBox
    SSTM_outsideMiddleCircle_downsamp = SSTM_outsideMiddleCircle(:,1:2:end);
    res15nm_outsideMiddleCircle = mean(SSTM_outsideMiddleCircle(:,[142 148]),2); 
    
    [recon_15pt_outsideMiddleCircle] = gaussSVD(SSTM_outsideMiddleCircle_downsamp,...
    res15nm_outsideMiddleCircle,sig);

% corner
    SSTM_corners_downsamp = SSTM_corners(:,1:2:end);
    res15nm_corners = mean(SSTM_corners(:,[142 148]),2); 
    
    [recon_15pt_corners] = gaussSVD(SSTM_corners_downsamp,...
    res15nm_corners,sig);
    
    
%% plot stuff

% Fig. 4(a)
figure;
subplot(2,2,1);
% imagesc(spectralRes_gaussFilt); colormap turbo
    
subplot(2,2,3)
    plot(spectralShift(1:20),SC_cropNorm_insideMiddleCircle);hold on
    plot(spectralShift(1:20),SC_cropNorm_outsideMiddleCircle);
    plot(spectralShift(1:20),SC_cropNorm_corners);
%     plot(spectralShift(1:20),SC_cropNorm_invertDonut);
    plot(spectralShift(1:20),ones(20,1)/2)
    axis([0 5 -.1 1])
    xlabel('Spectral Shift (nm)')
    ylabel('Norm. Spectral Correlation (C)')
    legend('Inner Circle','Outer Circle','Corners')

subplot(2,2,4)
    plot(calibrationWavelengths_fit(132:2:160,:),...
        recon_15pt_insideMiddleCircle(66:80));hold on
    plot(calibrationWavelengths_fit(132:2:160,:),...
        recon_15pt_outsideMiddleCircle(66:80));
    plot(calibrationWavelengths_fit(132:2:160,:),...
        recon_15pt_corners(66:80));
    axis([-inf inf -.05 .26])
    legend('Inner Circle','Outer Circle','Corners')
    xlabel('Wavelength (nm)')
    ylabel('Norm. Intensity')
    
    % visualize masks
figure;
    subplot(2,2,1)
    imagesc(mask_insideMiddleCircle);colormap gray
    subplot(2,2,2)
    imagesc(mask_outsideMiddleCircle);colormap gray
    subplot(2,2,3)
    imagesc(mask_corners);colormap gray

    
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


