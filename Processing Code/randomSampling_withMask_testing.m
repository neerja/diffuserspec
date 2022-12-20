% random sampling test with masking

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

%%

samplepercent = .1;
maskbox = 500;
[spectralPSFrand, samp_ind_nonan] = randSamplePSF_mask_jdm(spectralPSF_3D,samplepercent,maskbox);


spectralPSF_2D_downsamp = spectralPSF_2D(:,1:2:end);


    % grab 6 monochromatic frames evenly distributed over full bandwidth
    % note: each frame is an even indexed frame and not part of the downsampled
    % transfer matrix
    stepsize = ceil(344/6);
    offset = floor(round(344/6)-(round(344/6)/2));
    multiPeakVec = [stepsize-offset,(stepsize*2)-offset,(stepsize*3)-offset,...
        (stepsize*4)-offset,(stepsize*5)-offset,(stepsize*6)-offset];
    spectrumForRecon_multiWave = spectralPSF_2D(:,multiPeakVec);

    % set threshold for Truncated SVD recon and run recon and set max value of
    % recon to 1
    thresh = 10;
    [recon_tSVD] = TruncSVD(spectralPSF_2D_downsamp,...
        spectrumForRecon_multiWave,thresh);
    %recon_tSVD = abs(recon_tSVD/max(recon_tSVD(:)));

    [recon_gSVD] = gaussSVD(spectralPSF_2D_downsamp,...
        spectrumForRecon_multiWave,50);
    %recon_gSVD = abs(recon_gSVD/max(recon_gSVD(:)));


    figure;plot(calibrationWavelengths_fit(:,1:2:end)',recon_tSVD(1:172,:));
    title('truncated SVD - binary thresholding')
    xlabel('wavelegnth (nm)')

    figure;plot(calibrationWavelengths_fit(:,1:2:end)',recon_gSVD(1:172,:));
    title('truncated SVD - filtered thresholding')
    xlabel('wavelegnth (nm)')

    figure;plot(calibrationWavelengths_fit(:,62:85)',recon_tSVD(62:85,3));hold on
            plot(calibrationWavelengths_fit(:,62:85)',recon_gSVD(62:85,3));hold off

%%

%FIGURE 1b
    
    % select frames for 1nm and 2nm spacing (each wavelength step,
    % corresponds to .25 nm so 1nm step = 4 frame seperation, 2nm step = 8
    % frame seperation. Frames are averaged together to simulate
    % simultaneous acquisition
    res15nm = mean(spectralPSF_2D(:,[142 148]),2);
    
    % set threshold for Truncated SVD recon and run recon and set max value of
    % recon to 1    
    thresh = 25;
    [recon_2pt_tSVD] = TruncSVD(spectralPSF_2D_downsamp,...
        res15nm,thresh);
    recon_2pt_tSVD = abs(recon_2pt_tSVD/max(recon_2pt_tSVD(:)));

    sig = 500;
    [recon_2pt_gSVD] = gaussSVD(spectralPSF_2D_downsamp,...
        res15nm,sig);
    recon_2pt_gSVD = abs(recon_2pt_gSVD/max(recon_2pt_gSVD(:)));

figure;plot(calibrationWavelengths_fit(:,1:2:end)',recon_2pt_tSVD(1:172));hold on
            plot(calibrationWavelengths_fit(:,1:2:end)',recon_2pt_gSVD(1:172),'g');hold off

            %%

vertLines = zeros(size(spectralPSF_2D,2),6);
for vert = 1:6
    vertLines(multiPeakVec(vert)) = 1;
end

figure;
subplot(1,3,[1 2])
plot(calibrationWavelengths_fit(:,1:2:end)',recon_gSVD(1:172,:));hold on
stem(calibrationWavelengths_fit(:,multiPeakVec)',vertLines(vertLines == 1),'k--');hold off
axis([-inf inf 0 1.1])
xlabel('Wavelegnth (nm)')  
ylabel('Norm. Intensity')

subplot(1,3,3)
plot(calibrationWavelengths_fit(:,130:2:160)',recon_2pt_gSVD(65:80));hold on
stem(calibrationWavelengths_fit(:,[142 148])',[1 1],'k--');hold off
axis([818 825 0 1.1])
xlabel('Wavelegnth (nm)')  
ylabel('Norm. Intensity')