% Subfunctions: MakeDiffuserMatrix, simulateGaussianSource, generateNoise,
%               simulateOCTinterferogram

close all

%wave charicteristics
amplitude = 1;
nPixels = 1024;
lambda = 1.260E-6:.080E-6/(nPixels-1):1.340E-6; %(nm)
k = 2*pi./lambda; 
depth = [1 2];
z = [0 100e-6];

%initialize noise paramteters
stddev_laser = linspace(.0005,.05,1000);
stddev_detector = linspace(.0000005,.00005,1000);

%Make random diffuser matrix
diffuserMatrix = MakeDiffuserMatrix(nPixels);
    
%Initialize loop outputs
noisyInterferogram = zeros(numel(stddev_laser),nPixels);
spectrum_recon = zeros(nPixels,numel(stddev_laser));
    
for n = 1:numel(stddev_laser)
    
    %Simulate Gaussian Source
    [wavelength, intensityWavlength] = simulateGaussianSource(830E-9, 80E-9, nPixels);

    %Add source noise
    noisySource = generateNoise(intensityWavlength,0,stddev_laser(1));

    %Simulate OCT interferogram w/ n reflectors
    noisyInterferogram(n,:) = simulateOCTinterferogram(wavelength, noisySource, depth, z);
    
    %multiply diffuser and planewave of each wavelength
    detectedPSF = diffuserMatrix*noisyInterferogram(n,:)';

    %Add detector noise
    noisyPSF = generateNoise(detectedPSF',0,stddev_detector(n));

    %reconstruct spectrum using inverse diffuser matrix * detected PSF
    spectrum_recon(:,n) = diffuserMatrix\noisyPSF';

end
    
% Create A-scan from normal and reconstructed spectrum
Ascan = abs(ifft(noisyInterferogram'));
Ascan = Ascan(nPixels/2+1:end,:);

Ascan_recon = abs(ifft(spectrum_recon));
Ascan_recon = Ascan_recon(nPixels/2+1:end,:);

for scan = 1:n
    
    %Calculate SNR difference
    SNR_orig(scan) = 20*log10(max(Ascan(:,scan))/std(Ascan(300:400,scan)));
    SNR_recon(scan) = 20*log10(max(Ascan_recon(:,scan))/std(Ascan_recon(300:400,scan)));
    
end


%% Plots

figure;
plot(stddev_detector,SNR_orig); hold on; plot(stddev_detector,SNR_recon);hold off
title('SNR Comparison')
xlabel('Detector Noise (std)')
ylabel('SNR (dB)')
legend('SNR orig','SNR recon');

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

