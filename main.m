% main.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%
%  Purpose: 
%
%  Other m-files required: 

%  Subfunctions: 

%  MAT-files required: 

%--------------------------------------------------------------------------
%
%  Author:          Joe Malone (joseph.d.malone@vanderbilt.edu, Neerja
%                   Aggarwal (neerja@berkeley.edu)
%  Organization:    Bowden Biomedical Optics Lab, Vanderbilt University
%  Creation Date:   2022/20/12
%  Last Modified:   2022/20/12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%
% MAIN %
%%%%%%%%

%------------- BEGIN CODE --------------% 

%% 1. load calibration files

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

%% 2. load measurement

% ground truth spectrum collected on Thorlabs spectrometer 
load('./Data/broadbandMeasurements/Inphenix_groundTruth_wavelengths.mat');
load('./Data/broadbandMeasurements/Inphenix_groundTruth_spectrum.mat');

% assign filenames for measurement data
measurement_fileName = './Data/broadbandMeasurements/Inphenix_1.tif';
measurementBackground_fileName = './Data/broadbandMeasurements/broadband_openPinhole_background.tif';

% load diffuse broadband image of Inphenix source
broadbandMeasurement = tifRead(measurement_fileName);

% load background frame
broadbandMeasurement_background = tifRead(measurementBackground_fileName);

% perform background subtraction
broadbandMeasurement = abs(broadbandMeasurement - broadbandMeasurement_background);

%% 3. pre-process as needed

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
    b_broadbandIm = zeros(size(SSTM_outerRad,1),1);

    % Sample broadband diffuse image using coordinates from above
    for coordinate = 1:size(SSTM_outerRad,1)
        b_broadbandIm(coordinate,:) = broadbandIm(maskCoordinates_outerRad...
            (samp_corners(coordinate,1),2),maskCoordinates_outerRad...
            (samp_corners(coordinate,1),1),:);
    end



%% 4. run recon

    % assign sigma filter value and perform reconstruction
    sig =34;
    [recon_broadband] = gaussSVD(SSTM_outerRad(:,1:344),...
        b_broadbandIm,sig);

    % normalize to max value
    recon_broadband = recon_broadband/max(recon_broadband(10:330));


%% 5. analyze the spectral correlation for this diffuser matrix