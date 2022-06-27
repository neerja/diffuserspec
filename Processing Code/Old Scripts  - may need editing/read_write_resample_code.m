%load transfermatrix
transferMatrix_3D = uint16(tifRead_sequence('transferMatrix_3D',1,0));

tifWrite(transferMatrix_3D,'220grit_transferMatrix_double')


broadbandData_1 = uint16(mean(tifRead('broadband_diffuse_#0001.tif'),3));
tifWrite(broadbandData_1,'broadband_diffuse_#0001_uint16')
broadbandData_2 = uint16(mean(tifRead('broadband_diffuse_#0003.tif'),3));
tifWrite(broadbandData_2,'broadband_diffuse_#0002_uint16')
broadbandData_3 = uint16(mean(tifRead('broadband_diffuse_#0003.tif'),3));
tifWrite(broadbandData_3,'broadband_diffuse_#0003_uint16')


%% resample wavelength data

%fit calibration data with a 2nd order polynomial and plot
x = linspace(1,344,344);
[fit,pval] = polyfit(x,calibrationWavelengths,2);
calibrationWavelengths_polyfit = polyval(fit,x);

figure;plot(calibrationWavelengths);hold on;plot(calibrationWavelengths_polyfit);hold off
legend('Raw calib wavelengths','polynomial fit')
xlabel('pixel number')
ylabel('wavelength')

figure;plot(calibrationWavelengths_polyfit-calibrationWavelengths')
title('residual of polynomial fit')

%resample non-linear wavelength vector to a linear distribution
wavelength_vector = linspace(calibrationWavelengths(1),calibrationWavelengths(end),344);
calibrationWavelengths_resampled = interp1(calibrationWavelengths_polyfit...
    ,calibrationWavelengths_polyfit,wavelength_vector');

figure;plot(calibrationWavelengths_resampled);hold on;
    plot(calibrationWavelengths_polyfit);hold off
    legend('Resampled wavelengths','polynomial fit calib wavelenghts')
xlabel('pixel number')
ylabel('wavelength')

