%% spectralRecon_fig2v1.m
% Neerja Aggarwal & Joseph Malone
% July 13th, 2022
% Purpose: generate the recon for measurements.  

%% load calibration data

% load in stf matrix (randomly sampled):

% TODO:
% make a function that creates a sampling mask. takes in three parameters:
% example imagesize, min square,max square, number of sampling points

% for now, use the randomly sampled across the image.  
load 'Datasets matFiles/120grit/spectralPSF_randsamp_2D.mat';
% contains samp_line which is the actual PSFs
% samp_xy are the sampling indices to use for the test measurement.

% load wavelengths
load 'Datasets matFiles/calibrationFiles/calibrationWavelengths_fit.mat'

