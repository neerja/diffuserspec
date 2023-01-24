% Purpose: calculate the spectral correlation function based on equation 
% defined in DiffuserSpec manuscript, which was originally obtained from 
% [B. Redding, et al., Optics Express, 21, 5, 6584 (2013)] and 
% [P.Wang, R. Menon., Optics Express, 22,12,14575 (2014) 
% for a given spatial-spectral transfer matrix (SSTM) collected on the
% DiffuserSpec system.
% Inputs: 
%       SSTM = 2D NxM (spatialPixels x wavelength) matrix
%       iminWaveShift = starting index of number of wavelengths to evaluate. Note:
%           this is in pixels. For example, if the SSTM is 2000x100, where
%           100 wavelength pixels spans 25nm (1 pix. = 0.25nm), a maxWaveShift
%           of 20 will evaluate the spectral correlation over 20 pixels
%           which, corresponds to 20(0.25) = 5nm. 
%       imaxWaveShift = ending index of number of wavelengths to evaluate. Note:
%           this is in pixels. For example, if the SSTM is 2000x100, where
%           100 wavelength pixels spans 25nm (1 pix. = 0.25nm), a maxWaveShift
%           of 20 will evaluate the spectral correlation over 20 pixels
%           which, corresponds to 20(0.25) = 5nm. 
% Outpus: 
%   spectralCorrelation = 2D matrix containing computed correlation values
%   as a function of pixel and wavelengths. 

function [spectralCorrelation] = calcSpectralCorrelation(SSTM, iminWaveShift, imaxWaveShift)   

% truncate SSTM matrix to only include relevant wavelengths
SSTM_trunc = SSTM(:,iminWaveShift:imaxWaveShift);

% initialize SC variable
spectralCorrelation = zeros(size(SSTM_trunc,1),size(SSTM_trunc,2));

    for shift = 0:size(SSTM_trunc,2)-1
    
        % calculate SC based on equation defined in DiffuserSpec manuscript,
        % which was originally obtained from [B. Redding, et al., Optics
        % Express, 21, 5, 6584 (2013)] and [P.Wang, R. Menon., Optics Express,
        % 22,12,14575 (2014)
    
        s_row_trunc = SSTM_trunc(:,1:end-shift);
        s_row_shifted = circshift(SSTM_trunc,-shift,2);
        s_row_shifted_trunc = s_row_shifted(:,1:end-shift);
        
        SCtop = mean(s_row_trunc.*s_row_shifted_trunc,2);
        
        SCbottom = mean(s_row_trunc,2).*mean(s_row_shifted_trunc,2);
    
        spectralCorrelation(:,shift+1)= (SCtop./SCbottom) - 1;
    
    end

end