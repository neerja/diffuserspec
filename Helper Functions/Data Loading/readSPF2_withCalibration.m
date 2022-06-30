function [spectrum,wavelength] = readSPF2_withCalibration(calibrationCurve,...
    calibrationOrig,signalRange,filename)

    %if no filename is supplied, open the GUI to select the file
    if isempty(filename)
        disp('Select ground truth for measurement')
        [spf2File,path] = uigetfile('.spf2');
        filename = [path,spf2File]
    end
    
    data = readSPF2(filename);
    
    % if signalRange isn't specified, use default
    if isempty(signalRange)
        signalRange = [1740:2320];
    end
    
    wavelength = data(signalRange,1); %get subset of spectrometer input
    
    %apply correction spectrometer output and sample at calibration
    %wavelengths?  %ASK JOE
    spectrum = interp1(calibrationCurve,data(signalRange,2),calibrationOrig);
end