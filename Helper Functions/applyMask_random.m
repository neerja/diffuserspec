function [SSTM_2D,samp_ind] = applyMask_random(SSTM,maskCoordinates,samplepercent)
% This function's purpose is to randomly sample the spatio-spectral
% transfer matrix based on input percentage. 


% Inputs: 
%       SSTM = 2D NxM (spatialPixels x wavelength) matrix
%       maskCoordinates = coordinates of "on pixels" after mask is applied to SSTM.
%       samplepercent = scalar between 0 and 100.

    %initialize size of coordinate mask
    numCoords = size(maskCoordinates,1);
    
    %round the number of samples using the assinged percentage
    numsamples = floor(numCoords*(samplepercent/100));
    
    % create random pixel numbers within mask region to sample
    samp_ind = randi(numCoords,[numsamples,1]);

    % initialize variable
    SSTM_2D = zeros(numsamples,size(SSTM,3));

    % extract pixels from coordinate locations and concatonate into new
    % variable
    for coordinate = 1:numsamples
        SSTM_2D(coordinate,:) = SSTM(maskCoordinates(samp_ind(coordinate,1),2)...
            ,maskCoordinates(samp_ind(coordinate,1),1),:);

        % display progress
        if mod(coordinate,5000) == 0
            disp([num2str(coordinate)])
        else
        end

    end
end