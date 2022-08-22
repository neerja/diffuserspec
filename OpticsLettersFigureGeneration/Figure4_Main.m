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
    SSTM_2D = circshift(squeeze(SSTM(:,row,:)),172,2);

    for shift = 1:50

        % calculate spectral correlation numerator and denomenator (see
        % Menon paper for full expression)
        psftop = mean((SSTM_2D.*circshift(SSTM_2D,shift,2)),2);
        psfbottom = mean(SSTM_2D,2).*mean(circshift(SSTM_2D,shift,2),2);
    
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
%%%%% Inside of middle Box region %%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % assign percentage of data to sample and size of mask 
    samplepercent = 1;
    maskbox = 500;
    
    % perform random sampling over the specified mask region, region is
    % defined by the last input variable of randSamplePSF_mask_jdm
    % function, here 4 = inside of middle box mask
    [SSTM_insideMiddleBox, samp_insideMiddleBox,mask_insideMiddleBox]...
        = randSamplePSF_mask_jdm(SSTM,samplepercent,maskbox,4);
    
    % shift spectral dimension so that middle wavelength is first
    SSTM_insideMiddleBox = circshift(SSTM_insideMiddleBox,172,2);

    % calculate spectral correlation, limited to 50 wavelength shifts
    spectralCorrelation_insideMiddleBox = ...
        calcSpectralCorrelation(SSTM_insideMiddleBox, 50); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Outside of middle Box region %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    samplepercent = 0.22;
    maskbox = 500;

    [SSTM_outsideMiddleBox, samp_outsideMiddleBox,mask_outsideMiddleBox]...
        = randSamplePSF_mask_jdm(SSTM,samplepercent,maskbox,1);

    % shift spectral dimension so that middle wavelength is first
    SSTM_outsideMiddleBox = circshift(SSTM_outsideMiddleBox,172,2);

    spectralCorrelation_outsideMiddleBox = ...
        calcSpectralCorrelation(SSTM_outsideMiddleBox, 50);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%            Corners           %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    samplepercent = 1;
    maskbox = 500;

    [SSTM_corners, samp_corners,mask_corners] = ...
        randSamplePSF_mask_jdm(SSTM,samplepercent,maskbox,2);   

    % shift spectral dimension so that middle wavelength is first
    SSTM_corners = circshift(SSTM_corners,172,2);

    spectralCorrelation_corners = ...
        calcSpectralCorrelation(SSTM_corners, 50);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%         Invert Donut         %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    samplepercent = 0.22;
    maskbox = 500;
    [SSTM_invertdonut_singleout, samp_invertdonut_singleout,mask_invertdonut_singleout] = randSamplePSF_mask_jdm...
        (SSTM,samplepercent,maskbox,1);

    samplepercent = 4;
    maskbox = 250;
    [SSTM_invertdonut_doubleIn, samp_invertdonut_doubleIn,mask_invertdonut_doubleIn] = randSamplePSF_mask_jdm...
        (SSTM,samplepercent,maskbox,4);

    % Combine two mask regions
    samp_invertdonut = cat(1,samp_invertdonut_singleout,samp_invertdonut_doubleIn);
    SSTM_invertdonut = cat(1,SSTM_invertdonut_singleout,SSTM_invertdonut_doubleIn);
    
    % shift spectral dimension so that middle wavelength is first
    SSTM_invertdonut = circshift(SSTM_invertdonut,172,2);

    spectralCorrelation_invertdonut = ...
        calcSpectralCorrelation(SSTM_invertdonut, 50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Average along pixel locations and normalize each correlation function to 1
    SC_crop_insideMiddleBox = mean(spectralCorrelation_insideMiddleBox(:,1:20),1);
    SC_cropNorm_insideMiddleBox = SC_crop_insideMiddleBox/max(SC_crop_insideMiddleBox);

    SC_crop_outsideMiddleBox = mean(spectralCorrelation_outsideMiddleBox(:,1:20),1);
    SC_cropNorm_outsideMiddleBox = SC_crop_outsideMiddleBox/max(SC_crop_outsideMiddleBox);

    SC_crop_corners = mean(spectralCorrelation_corners(:,1:20),1);
    SC_cropNorm_corners = SC_crop_corners/max(SC_crop_corners);

    SC_crop_invertdonut = mean(spectralCorrelation_invertdonut(:,1:20),1);
    SC_cropNorm_invertDonut = SC_crop_invertdonut/max(SC_crop_invertdonut);

    
%%  Twopoint recons

%select sigma value for recon algoithm
    sig = 250;
    
% inside middleBox
    %downsample SSTM and select two wavelength frames for recon 
    SSTM__insideMiddleBox_downsamp = SSTM_insideMiddleBox(:,1:2:end);
    b_insideMiddleBox = mean(SSTM_insideMiddleBox(:,[142 148]),2); 

    [recon_15pt_insideMiddleBox] = gaussSVD(SSTM__insideMiddleBox_downsamp,...
    b_insideMiddleBox,sig);
    
% outsidemiddleBox
    SSTM_outsideMiddleBox_downsamp = SSTM_outsideMiddleBox(:,1:2:end);
    res15nm_outsideMiddleBox = mean(SSTM_outsideMiddleBox(:,[142 148]),2); 
    
    [recon_15pt_outsideMiddleBox] = gaussSVD(SSTM_outsideMiddleBox_downsamp,...
    res15nm_outsideMiddleBox,sig);

% corner
    SSTM_corners_downsamp = SSTM_corners(:,1:2:end);
    res15nm_corners = mean(SSTM_corners(:,[142 148]),2); 
    
    [recon_15pt_corners] = gaussSVD(SSTM_corners_downsamp,...
    res15nm_corners,sig);

% invert donut
    SSTM_invertdonut_conCat_double_downsamp = SSTM_invertdonut(:,1:2:end);
    res15nm_invertdonut_conCat_double = mean(SSTM_invertdonut(:,[142 148]),2); 
    
    [recon_15pt_invertdonut_conCat_double] = gaussSVD(SSTM_invertdonut_conCat_double_downsamp,...
    res15nm_invertdonut_conCat_double,sig);
    
    
%% plot stuff

% Fig. 4(a)
figure
imagesc(spectralRes_gaussFilt); colormap turbo
    
figure;
    plot(spectralShift(1:20),SC_cropNorm_insideMiddleBox);hold on
    plot(spectralShift(1:20),SC_cropNorm_outsideMiddleBox);
    plot(spectralShift(1:20),SC_cropNorm_corners);
    plot(spectralShift(1:20),SC_cropNorm_invertDonut);
    plot(spectralShift(1:20),ones(20,1)/2)
    axis([0 5 0 1])
    xlabel('Spectral Shift (nm)')
    ylabel('Norm. Spectral Correlation (C)')
    legend('Inner box','Outer box','Corners','Donut')

figure;
    plot(calibrationWavelengths_fit(1:2:end,:),recon_15pt_insideMiddleBox(1:end-1,:));hold on
    plot(calibrationWavelengths_fit(1:2:end,:),recon_15pt_outsideMiddleBox(1:end-1,:))
    plot(calibrationWavelengths_fit(1:2:end,:),recon_15pt_corners(1:end-1,:))
    plot(calibrationWavelengths_fit(1:2:end,:),recon_15pt_invertdonut_conCat_double(1:end-1,:))
    legend('inner box','outer box','corners','donut')
    xlabel('Wavelength (nm)')
    ylabel('Norm. Intensity')

    figure;
    plot(calibrationWavelengths_fit(132:2:160,:),...
        abs(recon_15pt_insideMiddleBox(66:80)));hold on
    plot(calibrationWavelengths_fit(132:2:160,:),...
        abs(recon_15pt_outsideMiddleBox(66:80)));
    plot(calibrationWavelengths_fit(132:2:160,:),...
        abs(recon_15pt_corners(66:80)));
    plot(calibrationWavelengths_fit(132:2:160,:),...
        abs(recon_15pt_invertdonut_conCat_double(66:80)));
    legend('Inner box','Outer box','Corners','Donut')
    xlabel('Wavelength (nm)')
    ylabel('Norm. Intensity')
    
    % visualize masks
figure;
    subplot(2,2,1)
    imagesc(mask_insideMiddleBox);colormap gray
    subplot(2,2,2)
    imagesc(mask_outsideMiddleBox);colormap gray
    subplot(2,2,3)
    imagesc(mask_corners);colormap gray
    subplot(2,2,4)
    imagesc(mask_invertdonut_singleout+mask_invertdonut_doubleIn);colormap gray

    legend('Inner box','Outer box','Corners','Donut')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spectralPSF_3D] = load_normalize_SSTM(spectralPSF_3D,backgroundFrame)

%subtract background and resign negative values
    spectralPSF_3D = abs(spectralPSF_3D - backgroundFrame);

% get average intensity for each spectral PSF and normalize max to 1
    normalizedIntensity = mean(mean(spectralPSF_3D,1),2);
    normalizedIntensity = normalizedIntensity/max(normalizedIntensity);

% normalize spectral intensity and set max value to 1
    spectralPSF_3D = spectralPSF_3D./normalizedIntensity;
    spectralPSF_3D = spectralPSF_3D/max(spectralPSF_3D(:));

end

function [spectralPSFrand, samp_ind_nonan,mask] = randSamplePSF_mask_jdm(spectralPSF_3D,samplepercent,maskbox,pattern)

% samples an image randomly.  sets values inside center box of size
% maskbox (half width) to 0 and omits those elements in output.

tic

    [N1,N2,N3] = size(spectralPSF_3D);

    numsamples = floor(N1*N2*(samplepercent/100));

    samp_ind = randi(N1*N2,[numsamples,1]);

    

    if pattern == 1 %middle omit

        mask = ones(N1,N2);
        % 2. get indices for mask
        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        mask(rows,cols) = 0;

    elseif pattern == 2 %corners

        mask = zeros(N1,N2);

        rows = [1:maskbox,N1-maskbox:N1];
        cols = [1:maskbox,N2-maskbox:N2];

        mask(rows,cols) = 1;

    elseif pattern == 3 % donut
        
        mask = zeros(N1,N2);

        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        rowsub = [floor(N1/2)-(maskbox/2):floor(N1/2)+(maskbox/2)];
        colsub = [floor(N2/2)-(maskbox/2):floor(N2/2)+(maskbox/2)];

        mask(rows,cols) = 1;
        mask(rowsub,colsub) = 0;

    elseif pattern == 4 %middle keep

        mask = zeros(N1,N2);
        % 2. get indices for mask
        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        mask(rows,cols) = 1;

    elseif pattern == 5 % invert donut
        
        mask = ones(N1,N2);

        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        rowsub = [floor(N1/2)-(maskbox/2):floor(N1/2)+(maskbox/2)];
        colsub = [floor(N2/2)-(maskbox/2):floor(N2/2)+(maskbox/2)];

        mask(rows,cols) = 0;
        mask(rowsub,colsub) = 1;

    elseif pattern == 6 % invert donut
        
        mask = ones(N1,N2);

        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        rowsub = [floor(N1/2)-(maskbox/4):floor(N1/2)+(maskbox/4)];
        colsub = [floor(N2/2)-(maskbox/4):floor(N2/2)+(maskbox/4)];

        mask(rows,cols) = 0;
        mask(rowsub,colsub) = 1;
    end

    % 3. set value to 0 inside A for mask
    spectralPSF_3D = spectralPSF_3D.*mask;

toc
    spectralPSF_3D = squeeze(reshape(spectralPSF_3D,N1*N2,1,N3));

    % check to see which ones are 0
    samp_ind_nonan = samp_ind(spectralPSF_3D(samp_ind) ~= 0);

    % take indexed nonNaN values for 2D PSF
    spectralPSFrand = spectralPSF_3D(samp_ind_nonan,:);
toc
end

function [spectralCorrelation] = calcSpectralCorrelation(SSTM, maxWaveShift)   
% This function's purpose is to calculate the spectral correlation function
% for a given spatial-spectral transfer matrix (SSTM) collected on the
% DiffuserSpec system.

% Inputs: 
%       SSTM = 2D NxM (spatialPixels x wavelength) matrix
%       maxWaveShift = index of number of wavelengths to evaluate. Note:
%           this is in pixels. For example, if the SSTM is 2000x100, where
%           100 wavelength pixels spans 25nm (1 pix. = 0.25nm), a maxWaveShift
%           of 20 will evaluate the spectral correlation over 20 pixels
%           which, corresponds to 20(0.25) = 5nm. 


    % initialize SC variable
    spectralCorrelation = zeros(size(SSTM,1),maxWaveShift);
    
    for shift = 1:maxWaveShift
    
        % calculate SC based on equation defined in DiffuserSpec manuscript,
        % which was originally obtained from [B. Redding, et al., Optics
        % Express, 21, 5, 6584 (2013)] and [P.Wang, R. Menon., Optics Express,
        % 22,12,14575 (2014)
        
        SCtop = mean((SSTM.*circshift(SSTM,shift,2)),2);
        SCbottom = mean(SSTM,2).*mean(circshift(SSTM,shift,2),2);
    
        spectralCorrelation(:,shift)= (SCtop./SCbottom) - 1;
    
    end
end

function [recon,dvec_inv,dvec_invf,gaussian_filter] = gaussSVD(A,b,sigma)

    %Inputs: A - 2D SSTM, obtained from calibration
    %        b - input diffuse spectrum to be reconstructed
    %        sigma - gaussian filter kernal for thresholding
    
    %Outputs: x = reconstructed input spectrum
    %       dvec_inv = inverse singular values
    %       dvec_invf = filtered inverse singular values
    %       gaussian_filter = soft threshold filter for inverse singular values
    
    % initialize size of input matrix A
        [numPixels,numWavelengths] = size(A);
    
    % perform SVD on SSTM
        [U,D,V] = svd(A);
    
    % create singular value vector by obtaining diagnol values of D
        dvec = diag(D);
        indexvec = 1:length(dvec);
    
    mu = 1; % (center for gauss filter)
    a = 1; % scaling coefficient
    
    % calculate normalized gaussian filter
        gaussian_filter = exp(-(indexvec-mu).^2/(2*sigma^2)) / (sigma*sqrt(2*pi));
        gaussian_filter = a*gaussian_filter/max(gaussian_filter); %normalize and scale
    
    % invert singular value vector and multiply by gaussian filter
        dvec_inv = 1./dvec;
        dvec_invf = dvec_inv.*gaussian_filter';
    
    % make diagnol matrix from inverse of D
        Df_inv = diag(dvec_invf);
        Z = zeros(numWavelengths,abs(numPixels-numWavelengths));
        Df_inv = [Df_inv Z]; 
    
    % solve for inverse of A matrix
        Ainvf =V*Df_inv*U';
        
    % perform reconstruction  
        recon = Ainvf*b;


end

