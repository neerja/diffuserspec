function [spectralCorrelation] = calcSpectralCorrelation(SSTM, maxWaveShift)   
% This function's purpose is to calculate the spectral correlation function
% for a given spatial-spectral transfer matrix (SSTM) collected on the
% DiffuserSpec system.

% Inputs: 
%       SSTM = 2D NxM (spatialPixels x wavelength) matrix
%       maxWaveShift = index of number of wavelengths to evaluate. Note:
%           this is in pixels. For example, if the SSTM is 2000x100, where
%           100 wavelength pixels spans 25nm (1 pix. = 0.25nm), a maxWaveShift
%           of 20 will evaluate the spectral correlation over 20 pixels
%           which, corresponds to 20(0.25) = 5nm. 


% initialize SC variable
spectralCorrelation = zeros(size(SSTM,1),maxWaveShift);

    for shift = 1:maxWaveShift
    
        % calculate SC based on equation defined in DiffuserSpec manuscript,
        % which was originally obtained from [B. Redding, et al., Optics
        % Express, 21, 5, 6584 (2013)] and [P.Wang, R. Menon., Optics Express,
        % 22,12,14575 (2014)
    
        SCtop = mean((SSTM.*circshift(SSTM,shift,2)),2);
        SCbottom = mean(SSTM,2).*mean(circshift(SSTM,shift,2),2);
    
        spectralCorrelation(:,shift)= (SCtop./SCbottom) - 1;
    
    end

end