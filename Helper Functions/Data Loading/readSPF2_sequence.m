function [spectrum,wavelength] = readSPF2_sequence(calibrationCurve,...
    calibrationOrig,signalRange,varargin)

% if no input - code will ask you to select folder
% otherwise, varargin = folder name you want to process

oldfolder = cd;

if isempty(varargin)
    folderPath = uigetdir;
else
    folderPath = varargin{1};
end
cd(folderPath)

spf2Files = dir('*.spf2');
temp = readSPF2(spf2Files(1).name);
wavelength = temp(signalRange,1);

spectrum = zeros(size(temp(signalRange,2),1),size(spf2Files,1));
    for n = 1:numel(spf2Files)

        temp = readSPF2(spf2Files(n).name);
        %apply calibration data to spectrometer output
        spectrum(:,n) = interp1(calibrationCurve,temp(signalRange,2),calibrationOrig);
    
    end
    
    cd(oldfolder)
    
end