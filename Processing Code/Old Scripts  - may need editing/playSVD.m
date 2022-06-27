%% playSVD.m
% Neerja Aggarwal
% play with truncated SVD.  Why is there an oscillating behavior in
% truncated SVD? 

% import calibration matrix 
%% load diffuser spec data

%load calibration matrix
disp('Select calibration matrix file')
spectrumForRecon = uigetfile('.tif');
spectralPSF_3D = tifRead(spectrumForRecon);

%load background frame
disp('Select background file')
backgroundFilename = uigetfile('.tif');
backgroundFrame = tifRead(backgroundFilename);

%subtract background and resign negative values
spectralPSF_3D = abs(spectralPSF_3D - backgroundFrame);
% get rid of last frame cause it's empty
spectralPSF_3D = spectralPSF_3D(:,:,1:end-1);

%% Spectral normalization of calibration data

% get average intensity for each spectral PSF and set range between [0,1]
normalizedIntensity = mean(mean(spectralPSF_3D,1),2);
normalizedIntensity = normalizedIntensity/max(normalizedIntensity);

% normalize spectral intensity and set max value to 1
spectralPSF_3D = spectralPSF_3D./normalizedIntensity;
spectralPSF_3D = spectralPSF_3D/max(spectralPSF_3D(:));

%% make 2D

line_number = 900;
A = squeeze(spectralPSF_3D(line_number,:,:));

%% A = calibrated diffuser transform matrix
[m,n] = size(A);

[U,D,V] = svd(A);
%%
dvec = diag(D);
figure
semilogy(dvec)
xlabel('Index')
ylabel('Singular Value')

dvecinv = 1./dvec;
figure
semilogy(dvecinv)
xlabel('Index')
ylabel('Inverse Singular Value')
%%

thresh = 2;
% D = diagnol matrix containing singular values
Dprime = (1./D)';
Dprime(Dprime == inf) = 0;
Dprime_thresh = Dprime;


% set any singular values (inverted) above threshold to zero
Dprime_thresh(Dprime_thresh > thresh) = 0;

%% make low rank A:
thresh = 0.1;
Dthresh = D;
Dthresh(Dthresh<1/thresh) = 0;
Alr = U*Dthresh*V';
imagesc(Alr)
%%
y = V'*recon;