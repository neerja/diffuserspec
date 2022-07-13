% processDiffuserSpec_MAIN.m

% This code is meant to serve as the main processing code for the
% diffuserSpec manuscript submitted to Optics Letters 

% this code will read in the spectral calibration data for a given diffuser
% and run various scripts to generate figures that can be found in the
% manuscript


% Author - Joe Malone (joseph.d.malone@vanderbilt.edu)
% Most Recent Edits - 06/29/2022

%% IMPORTANT - Starting Directory

% Note: the data required for this code is located on the Projects server 
% of the BBOL lab, depending on the drive location on your local PC, 
% you may need to appropriate drive. Also, the data files required for 
% processing are rather large, so it is recommended to save the files on
% your local drive

% In the case that you save the data from the server on your local drive,
% select the local folder as the originalDir, otherwise use the specified
% server location

% SET DIR to - Projects\smartOCT\Analysis\Data\CDMRPVRP\diffuserSpecData_OpticsLetters\organizedCodeData\Data
originalDir = uigetdir;

%% load in pre calibrated data and constants 

% select calibrationFiles folder (located with
calibrationDir = uigetdir;
cd(calibrationDir)

load('calibrationOrig.mat');
load('calibrationCurve.mat');
load('signalRange.mat');
load('calibrationWavelengths_fit.mat');
load('partialSpectrumWavelengths_gt');
load('partialSpectrum_broadband1');
load('partialSpectrum_broadband2');
load('partialSpectrum_gt_2');
load('partialSpectrum_largeBroadband.mat')

cd(originalDir)

%% load diffuser spec data

%load calibration matrix
disp('Loading 120grit diffuser calibration matrix')
spectralPSF_3D = tifRead('120grit_spectralPSF_3D.tif');

%load background frame
disp('Loading background signal for background subtraction')
backgroundFrame = tifRead('120grit_spectralPSF_background.tif');

%subtract background and resign negative values
spectralPSF_3D = abs(spectralPSF_3D - backgroundFrame);

% get average intensity for each spectral PSF and normalize max to 1
normalizedIntensity = mean(mean(spectralPSF_3D,1),2);
normalizedIntensity = normalizedIntensity/max(normalizedIntensity);

% normalize spectral intensity and set max value to 1
spectralPSF_3D = spectralPSF_3D./normalizedIntensity;
spectralPSF_3D = spectralPSF_3D/max(spectralPSF_3D(:));

clear backgroundFrame normalizedIntensity

%% Generate Figure 2

% FIGURE 1a

    % Create 2D spectral transfer function using line from the middle of the 2D
    % frame, line 1080
    lineNum = 1080;
    spectralPSF_2D_full_1080 = squeeze(spectralPSF_3D(lineNum,:,1:344));
    
    % downsample 2D transfer function by a factor of 2 for reconstruction, code
    % takes every other frame starting at 1 (all odd frames)
    downSampFac = 2;
    spectralPSF_2D_downSamp_1080 = spectralPSF_2D_full_1080(:,1:downSampFac:end);
    
    % grab 6 monochromatic frames evenly distributed over full bandwidth
    % note: each frame is an even indexed frame and not part of the downsampled
    % transfer matrix
    stepsize = ceil(344/6);
    offset = floor(round(344/6)-(round(344/6)/2));
    multiPeakVec = [stepsize-offset,(stepsize*2)-offset,(stepsize*3)-offset,...
        (stepsize*4)-offset,(stepsize*5)-offset,(stepsize*6)-offset];
    spectrumForRecon_multiWave = spectralPSF_2D_full_1080(:,multiPeakVec);
    
    % set threshold for Truncated SVD recon and run recon and set max value of
    % recon to 1
    thresh = 16;
    [recon] = TruncSVD(spectralPSF_2D_downSamp_1080,...
        spectrumForRecon_multiWave,thresh);
    recon_full = abs(recon/max(recon(:)));

%FIGURE 1b
    
    % select frames for 1nm and 2nm spacing (each wavelength step,
    % corresponds to .25 nm so 1nm step = 4 frame seperation, 2nm step = 8
    % frame seperation. Frames are averaged together to simulate
    % simultaneous acquisition
    res2nm = mean(spectralPSF_2D_full_1080(:,[172 180]),2);
    res15nm = mean(spectralPSF_2D_full_1080(:,[172 178]),2);
    
    % set threshold for Truncated SVD recon and run recon and set max value of
    % recon to 1
    thresh = 25;
    [recon] = TruncSVD(spectralPSF_2D_downSamp_1080,...
        res2nm,thresh);
    recon_2nm_1080 = abs(recon/max(recon(:)));
    
    thresh = 25;
    [recon] = TruncSVD(spectralPSF_2D_downSamp_1080,...
        res15nm,thresh);
    recon_15nm_1080 = abs(recon/max(recon(:)));

%FIGURE 1c

    % Create 2D spectral transfer function using different line of the 2D
    % frame, line 400
    lineNum = 400;
    spectralPSF_2D_full_400 = squeeze(spectralPSF_3D(lineNum,:,1:344));

    % downsample 2D transfer function by a factor of 2 for reconstruction, code
    % takes every other frame starting at 1 (all odd frames)
    spectralPSF_2D_downSamp_400 = spectralPSF_2D_full_400(:,1:downSampFac:end);

    res15nm_400 = mean(spectralPSF_2D_full_400(:,[172 178]),2);

    thresh = 25;
    [recon,~,Dprime] = TruncSVD(spectralPSF_2D_downSamp_400,...
        res15nm_400,thresh);
    recon_15nm_400 = abs(recon/max(recon(:)));



%% Plot Fig2 data
figure;
subplot(2,3,[1:3])
plot(calibrationWavelengths_fit(:,1:downSampFac:end)',recon_full);
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

subplot(2,3,5)
plot(calibrationWavelengths_fit(:,1:downSampFac:end)',recon_15nm_1080);
axis([822 838 0 1])
title('Row 1080 - 1.5 nm')
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

subplot(2,3,4)
plot(calibrationWavelengths_fit(:,1:downSampFac:end)',recon_2nm_1080);
axis([822 838 0 1])
title('Row 1080 - 2 nm')
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

subplot(2,3,6)
plot(calibrationWavelengths_fit(:,1:downSampFac:end)',recon_15nm_400);
axis([822 838 0 1])
title('Row 400 - 1.5 nm')
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')

%% 

% FIGURE 3 Broadband data


lineNum = 400;

% load broadband spectrum
spectrumForRecon1 = tifRead('120grit_broadband_1sthalf.tif');
spectrumForRecon2 = tifRead('120grit_broadband_2ndhalf.tif');
spectrumForRecon1_400 = spectrumForRecon1(lineNum,:)';
spectrumForRecon2_400 = spectrumForRecon2(lineNum,:)';

spectrumForReconFull = (spectrumForRecon1+spectrumForRecon2)/2;
spectrumForRecon_400 = spectrumForReconFull(lineNum,:)';

spectrumForRecon_partial = tifRead('120grit_broadband_#0001.tif');
spectrumForRecon_partial_400 = spectrumForRecon_partial(lineNum,:)';


thresh = 1;
[recon,~,Dprime,Dprime_thresh] = TruncSVD(spectralPSF_2D_full_400,...
    spectrumForRecon1_400,thresh);
recon_Full1 = abs(recon/max(recon(:)))-min(abs(recon/max(recon(:))));

thresh = 1;
[recon,~,Dprime,Dprime_thresh] = TruncSVD(spectralPSF_2D_full_400,...
    spectrumForRecon2_400,thresh);
recon_Full2 = abs(recon/max(recon(:)))-min(abs(recon/max(recon(:))));


thresh = 1.2;
[recon,~,Dprime,Dprime_thresh] = TruncSVD(spectralPSF_2D_full_400,...
    spectrumForRecon_partial_400,thresh);
recon_partial = abs(recon/max(recon(:)))-min(abs(recon/max(recon(:))));

figure;
subplot(1,2,1);
plot(calibrationWavelengths_fit',recon_Full1(1:344,:));hold on
plot(calibrationWavelengths_fit',recon_Full2(1:344,:));
plot(partialSpectrumWavelengths_gt(28:541),partialSpectrum_broadband1(28:541))
plot(partialSpectrumWavelengths_gt(28:541),partialSpectrum_broadband2(28:541))
axis([-inf inf 0 1])
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')
legend('Recon.','GT')

subplot(1,2,2)
plot(calibrationWavelengths_fit',recon_partial(1:344,:));hold on
plot(partialSpectrumWavelengths_gt,partialSpectrum_largeBroadband)
axis([-inf inf 0 1])
xlabel('Wavelength (nm)')
ylabel('Norm. Intensity')
legend('Recon.','GT')


%% Figure 3

refWave = [25,172,319];
refWave_realWaves = [calibrationWavelengths_fit(15),...
    calibrationWavelengths_fit(172),calibrationWavelengths_fit(329)];

spectralPSF_2D_full_400 = spectralPSF_2D_full_400';
spectralCorrelation_Normalized = zeros(344,6);

    % spectral correlation vs wavelength number 1
    PSFt = mean(spectralPSF_2D_full_400(refWave(1),:).*spectralPSF_2D_full_400,2);
    PSFb = mean(spectralPSF_2D_full_400(refWave(1),:),2).*mean(spectralPSF_2D_full_400,2);
    spectralCorrelation = (PSFt./PSFb)-1;
    spectralCorrelation_Normalized(:,1) = ...
        spectralCorrelation/max(spectralCorrelation);

    % spectral correlation vs wavelength number 2
    PSFt = mean(spectralPSF_2D_full_400(refWave(2),:).*spectralPSF_2D_full_400,2);
    PSFb = mean(spectralPSF_2D_full_400(refWave(2),:),2).*mean(spectralPSF_2D_full_400,2);
    spectralCorrelation = (PSFt./PSFb)-1;
    spectralCorrelation_Normalized(:,2) = ...
        spectralCorrelation/max(spectralCorrelation);

    % spectral correlation vs wavelength number 3
    PSFt = mean(spectralPSF_2D_full_400(refWave(3),:).*spectralPSF_2D_full_400,2);
    PSFb = mean(spectralPSF_2D_full_400(refWave(3),:),2).*mean(spectralPSF_2D_full_400,2);
    spectralCorrelation = (PSFt./PSFb)-1;
    spectralCorrelation_Normalized(:,3) = ...
        spectralCorrelation/max(spectralCorrelation);

spectralPSF_2D_full_1080 = spectralPSF_2D_full_1080';

    % spectral correlation vs wavelength number 1
    PSFt = mean(spectralPSF_2D_full_1080(refWave(1),:).*spectralPSF_2D_full_1080,2);
    PSFb = mean(spectralPSF_2D_full_1080(refWave(1),:),2).*mean(spectralPSF_2D_full_1080,2);
    spectralCorrelation = (PSFt./PSFb)-1;
    spectralCorrelation_Normalized(:,4) = ...
        spectralCorrelation/max(spectralCorrelation);

    % spectral correlation vs wavelength number 2
    PSFt = mean(spectralPSF_2D_full_1080(refWave(2),:).*spectralPSF_2D_full_1080,2);
    PSFb = mean(spectralPSF_2D_full_1080(refWave(2),:),2).*mean(spectralPSF_2D_full_1080,2);
    spectralCorrelation = (PSFt./PSFb)-1;
    spectralCorrelation_Normalized(:,5) = ...
        spectralCorrelation/max(spectralCorrelation);

    % spectral correlation vs wavelength number 3
    PSFt = mean(spectralPSF_2D_full_1080(refWave(3),:).*spectralPSF_2D_full_1080,2);
    PSFb = mean(spectralPSF_2D_full_1080(refWave(3),:),2).*mean(spectralPSF_2D_full_1080,2);
    spectralCorrelation = (PSFt./PSFb)-1;
    spectralCorrelation_Normalized(:,6) = ...
        spectralCorrelation/max(spectralCorrelation);


    %% calculate spectral resolution vs line number 

    spectralCorrelation_Normalized_lineNum = zeros(344,size(spectralPSF_3D,1));
for lineNum = 1:size(spectralPSF_3D,1)

    spectralPSF_2D = squeeze(spectralPSF_3D(lineNum,:,1:344))';

    PSFt = mean(spectralPSF_2D(refWave(2),:).*spectralPSF_2D,2);
    PSFb = mean(spectralPSF_2D(refWave(2),:),2).*mean(spectralPSF_2D,2);
    spectralCorrelation = (PSFt./PSFb)-1;

    spectralCorrelation_Normalized_lineNum(:,lineNum) = ...
        spectralCorrelation/max(spectralCorrelation);

    clear spectralPSF_2D PSFt PSFb spectralCorrelation 

    lineNum
end

spectralRes = zeros(size(spectralPSF_3D,1),1);
thresh = 0.5;
for lineNum = 1:size(spectralPSF_3D,1)
    [ind,] = find(spectralCorrelation_Normalized_lineNum(172:end,lineNum)<thresh);

    if isempty(ind)
        spectralRes(lineNum) = spectralShift(1,173);
    else
        spectralRes(lineNum) = spectralShift(1,ind(1));
    end
end


    %% plot Fig. 3

temp = spectralPSF_2D_full_400/max(spectralPSF_2D_full_400(:));
figure;
subplot(2,2,2)
imagesc(calibrationWavelengths_fit,1:2560,temp,[-.05 1])
colormap hot
colorbar('eastoutside')
xlabel('Wavelength(nm)')
ylabel('X-Pixel Value')

waveRange = calibrationWavelengths_fit(end)-calibrationWavelengths_fit(1);
spectralShift = linspace(0,waveRange,344);
subplot(2,2,3)
plot(spectralShift(1:25),...
    spectralCorrelation_Normalized(refWave(1):(refWave(1)+24),1)); hold on
plot(spectralShift(1:25),...
    spectralCorrelation_Normalized(refWave(2):(refWave(2)+24),2));
plot(spectralShift(1:25),...
    spectralCorrelation_Normalized(refWave(3):(refWave(3)+24),3)); 

plot(spectralShift(1:25),...
    spectralCorrelation_Normalized(refWave(1):(refWave(1)+24),4));
plot(spectralShift(1:25),...
    spectralCorrelation_Normalized(refWave(2):(refWave(2)+24),5));
plot(spectralShift(1:25),...
    spectralCorrelation_Normalized(refWave(3):(refWave(3)+24),6));hold off

legend(['line 400 - ',num2str(refWave_realWaves(1))],['line 400 - ',num2str(refWave_realWaves(2))],...
    ['line 400 - ',num2str(refWave_realWaves(3))],['line 1080 - ',num2str(refWave_realWaves(1))],['line 1080 - ',num2str(refWave_realWaves(2))],...
    ['line 1080 - ',num2str(refWave_realWaves(3))])
xlabel('Spectral Shift (nm)')
ylabel('Norm. Correlation')

pixRange = 1:2160;
subplot(2,2,4)
plot(pixRange,spectralRes)
xlabel('Row Number')
ylabel('Spectral Resolution (nm)')



