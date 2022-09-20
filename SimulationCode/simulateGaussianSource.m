function [wavelength, gaussianSource] = simulateGaussianSource(centerWavelength, fullBandwidth, numPoints)


%% Laser source parameters
FWHMWavlength = 30e-9;

%% Simulate a gaussian source
wavelength = linspace(centerWavelength - fullBandwidth/2, centerWavelength + fullBandwidth/2, numPoints);
sigmaWavlength = FWHMWavlength / (2 * sqrt(2 * log(2)));
gaussianSource = gaussmf(wavelength, [sigmaWavlength, centerWavelength]);

%Display spectrum
% figure();
% plot(wavelength * 1e9, gaussianSource);
% title('Light source spectrum in wavelength domain');
% xlabel('Wavelength(nm)');
% ylabel('Intensity(a.u.)');

end