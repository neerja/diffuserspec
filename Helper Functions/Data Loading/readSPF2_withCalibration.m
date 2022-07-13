% readSPF2_withCalibration.m
% last argument (wavelengthInterp) is optional. It will resample the
% spectrum after calibration correction. 

function [spectrum,wavelength] = readSPF2_withCalibration(wavelengthOrig,...
    wavelengthCorrected,signalRange,filename1,wavelengthInterp)

    %if no filename is supplied, open the GUI to select the file
    if nargin<4
        disp('Select ground truth for measurement')
        [spf2File,path] = uigetfile('.spf2');
        filename1 = [path,spf2File]
    end
    
    data = readSPF2(filename1);
    
    % if signalRange isn't specified, use default
    if isempty(signalRange)
        signalRange = [1740:2320];
    end
    
    wavelength = data(signalRange,1); %get subset of spectrometer input
    
    %apply correction spectrometer output and sample at calibration (2 nm
    %shift)
    spectrum = interp1(wavelengthOrig,data(signalRange,2),wavelengthCorrected);

    if nargin<5
        %do nothing
    else
        spectrum = interp1(wavelength,spectrum,wavelengthInterp);
        wavelength = wavelengthInterp;
    end
end