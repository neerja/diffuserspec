%NEED Fuzzy Logic Toolbox

% close all

%wave charicteristics
amplitude = 1;
nPixels = 1024;
lambda = 1.260E-6:.080E-6/(nPixels-1):1.340E-6; %(nm)
k = 2*pi./lambda; 
n = [1 1.5 2 1];
z = [0 100e-6 200e-6 300e-6];

%Simulate Gaussian Source
[wavelength, intensityWavlength] = simulateGaussianSource(830E-9, 80E-9, nPixels);

%Add source noise
[noisySource] = generateNoise(intensityWavlength,0,.05);

%Simulate OCT interferogram w/ n reflectors
[noisyInterferogram] = simulateOCTinterferogram(wavelength, noisySource, n, z);

Ascan = abs(fftshift(fft(noisyInterferogram)));
Ascan = Ascan(nPixels/2+1:end);

%create random diffuser matrix
diffuserMatrix = MakeDiffuserMatrix(nPixels);

%multiply diffuser and planewave of each wavelength
detectedPSF = diffuserMatrix*noisyInterferogram';

%Add detector noise
[noisyPSF] = generateNoise(detectedPSF',0,.0001);

%reconstruct spectrum using inverse diffuser matrix * detected PSF
spectrum_recon = diffuserMatrix\noisyPSF';

Ascan_recon = abs(fftshift(fft(spectrum_recon)));
Ascan_recon = Ascan_recon(nPixels/2+1:end)';

%Calculate SNR difference
SNR_orig = 20*log10(max(Ascan)/std(Ascan(:,300:400)));
SNR_recon = 20*log10(max(Ascan_recon)/std(Ascan_recon(:,300:400)));

%% Plots

figure();
subplot(2,3,1); plot(wavelength,noisySource);
title('Source Spectrum')
xlabel('wavelength')
ylabel('Int.')

subplot(2,3,2); plot(wavelength,noisyInterferogram);
title('Interferogram')
xlabel('wavelength')
ylabel('Int.')

subplot(2,3,3); plot(Ascan);
title('A-scan')
xlabel('pixel')
ylabel('Int.')

subplot(2,3,4); imagesc(diffuserMatrix);
title('Diffuser Matrix')
xlabel('wavelength')
ylabel('pixel')

subplot(2,3,5); plot(wavelength, spectrum_recon);
title('Reconstructed interferogram');
xlabel('wavelength')
ylabel('Int.')

subplot(2,3,6); plot(Ascan_recon);
title('Reconstructed A-scan');
xlabel('pixel')
ylabel('Int')

%% Gradient Decent Reconstruction

% [im_reconstruct] = gradDecent_reconstruct(detectedPSF, interferogram, 1, 200);







