function [spectrum,wavelength] = readSPF2_withCalibration(calibrationCurve,...
    calibrationOrig,signalRange)

    [spf2File,path] = uigetfile('.spf2');
    origPath = cd;
    cd(path)

    data = readSPF2(spf2File);
    
    wavelength = data(signalRange,1);
    
    %apply calibration data to spectrometer output
    spectrum = interp1(calibrationCurve,data(signalRange,2),calibrationOrig);

    cd(origPath)
end