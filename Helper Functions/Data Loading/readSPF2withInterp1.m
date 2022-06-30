function [wv, spectrum] = readSPF2withInterp1(wavelengths,filename)

    %if no filename is supplied, open the GUI to select the file
    if isempty(filename)
        disp('Select ground truth for measurement')
        [spf2File,path] = uigetfile('.spf2');
        filename = [path,spf2File]
    end
    
    data = readSPF2(filename);
    wv = data(:,1);
    spectrum = data(:,2);
        
    if isempty(wavelengths)
        % do nothing
    else
        % make sure to run interp before overriding the wv
        spectrum = interp1(wv,spectrum,wavelengths);
        wv = wavelengths;
    end  
    
end
