%% spectralRecon_fig2v1.m
% Neerja Aggarwal & Joseph Malone
% July 13th, 2022
% Purpose: generate the recon for measurements.  

%% load calibration dataset (takes a minute)
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
spectralPSF_3D = tifRead(filename);

infofilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat'
load('./Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/calibrationInfo.mat')
% contains calibrationWavelengths_fit, description

%% TODO:
% make a function that creates a sampling mask. takes in three parameters:
% example imagesize, min square,max square, number of sampling points. 
% use randSamplePSF then ignore any values inside min square. 

%% Create randomly sampled 2D spectralPSF
% choose percentage of data to sample
samplepercent = 1; %percent of pixels
[spectralPSF_2D,sampxy,sampfac] = randSamplePSF(spectralPSF_3D,samplepercent);

%% load midband spectrum
filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/partialSpectrum_#0001.tif';
bgfilename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/partialSpectrum/background_partialSpectrum_#0001.tif';

spectrumForRecon_sampled = loadTestMeasurement(filename,bgfilename,sampxy,sampfac);

% load wavelengths
load 'Datasets matFiles/calibrationFiles/calibrationWavelengths_fit.mat'

%%

