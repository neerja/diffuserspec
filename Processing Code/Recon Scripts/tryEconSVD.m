% tryEconSVD.m

% try to use all the pixels.  Use econ SVD to compute the svd (will atleast
% give you 344 singular values & vectors).


% CHANGE DIRECTORY TO TOP LEVEL /DIFFUSERSPEC/
% load spectralSSTM matrix as A.  A is mxn where m>n. 

filename = './Raw Data/120grit_partialBroadbandSpectrum_2022-06-16/spectralPSF_3D_norm.tif'
A = tifRead(filename);

%% (optional) reduce size of A by half to make it easier to work 
% I only had a 16Gb laptop. (A=8 Gb instead of 15 Gb so then I can compare
% A-Ap
A = A(:,:,1:2:end); %every other wavelength

%% make 2D (~5 million pixels by 172 wavelength)
[N1,N2,N3]= size(A)
A = reshape(A,[N1*N2,N3]);

%% compute econ svd (doesn't compute unnecessary U vectors)
[U,dvec,V] = svd(A,"econ",'vector');

%% make Aprime from computed singular values & vectors. 
D = diag(dvec);
Ap = U*D*V';

%% compute the fro norm.
An = norm(A,'fro')
Apn = norm(A-Ap,'fro')
% An =    5.1334e+05
% Apn =   1.6800e-09
% basically using the econ SVD gives back all the same info to reconstruct
% the same A matrix. 

%% run inverse SVD code (compute inverse)
[m,n] = size(A);
dvec_inv = 1./dvec;
Df_inv = diag(dvec_inv);
% do not pad Dinv
% optional - filter Df_inv however you wish (gauss or box or whatever).
% Here I don't filter for demo purposes. 
Ainv =V*Df_inv*U';
% We are able to compute Ainv!

