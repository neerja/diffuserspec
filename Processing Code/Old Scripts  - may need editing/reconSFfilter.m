%% reconSVfilter.m
% Neerja Aggarwal
% June 20th, 2022
% Purpose: load data, take SVD, apply filter to singular values, use
% inverse to do recon.


% A = calibrated diffuser transform matrix, m>n
[m,n] = size(A);
[U,D,V] = svd(A);
%%
dvec = diag(D);
figure
semilogy(dvec)
xlabel('Index')
ylabel('Singular Value')

%% make gauss filter

sigma = 30 ; % (std dev for gauss filter)
indexvec = 1:length(dvec);
mean = 1; % (attenuation for gauss filter)

gaussian_filter = exp(-indexvec.^2/(2*sigma^2)) / (sigma*sqrt(2*pi));
gaussian_filter = gaussian_filter/max(gaussian_filter);

figure
plot(gaussian_filter)
xlabel('Index')

dvec_inv = 1./dvec;
dvec_invf = dvec_inv.*gaussian_filter';
dvec_invt = dvec_inv;
dvec_invt(dvec_invt>2) = 0; %appl

figure
semilogy(dvec_inv)
hold on
semilogy(dvec_invt)
semilogy(dvec_invf)
xlabel('Index')
ylabel('Singular Value')
legend('Original Dinv', 'Thresholded Dinv','Filtered Dinv')

figure
plot(dvec_inv)
hold on
plot(dvec_invt)
plot(dvec_invf)
xlabel('Index')
ylabel('Singular Value')
legend('Original Dinv', 'Thresholded Dinv','Filtered Dinv')

% make inverse and do recon

Df_inv = diag(dvec_invf);
Z = zeros(n,abs(m-n));
Df_inv = [Df_inv Z]; 
% apply threshold
% Df_inv(Df_inv > 10 ) = 0;

Ainvf =V*Df_inv*U';

% b = spectrumForRecon_sampled; 
b = b3;
x = Ainvf*b;
x = x./max(x);

figure
plot(calibrationWavelengths_fit,x)
hold on;
xlabel('Wavelength (nm)')
ylabel('Intensity')
%
thresh = 2;
x = TruncSVD(A,b,thresh);
x = x./max(x);

% plot(calibrationWavelengths_fit,x)

%
disp('Select ground truth for measurement')
[partialSpectrum_gt,partialSpectrumWavelengths_gt] = ...
    readSPF2_withCalibration(calibrationCurve,calibrationOrig,signalRange);
partialSpectrum_gt = partialSpectrum_gt/max(partialSpectrum_gt(:));

partialSpectrum_gt = imresize(partialSpectrum_gt,...
    [344 size(partialSpectrum_gt,2)]);
partialSpectrumWavelengths_gt = imresize(partialSpectrumWavelengths_gt,...
    [344 size(partialSpectrumWavelengths_gt,2)]);

plot(partialSpectrumWavelengths_gt,partialSpectrum_gt)
xlabel('Wavelength (nm)')
ylabel('normallized intensity')
% legend('gaussian recon','threshold recon', 'groundtruth')
legend('gaussian recon','groundtruth')
