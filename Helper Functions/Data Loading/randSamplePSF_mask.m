function [spectralPSFrand, samp_ind_nonan] = randSamplePSF_mask(spectralPSF_3D,samplepercent,maskbox)

% samples an image randomly.  sets values inside center box of size
% maskbox (half width) to NaN and omits those elements in output.

% 1. get random samples over linear indexing. 
    if nargin<3
        maskbox = 500; %half width of the box. 
    end
    
    if nargin<2
        samplepercent = 0.1;
    end
    
    [N1,N2,N3] = size(spectralPSF_3D);
    numsamples = floor(N1*N2*(samplepercent/100));
    samp_ind = randi(N1*N2,[numsamples,1]);
    
    % 2. get indices for mask
    rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
    cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];
    
    % 3. set value to NaN inside A for mask
    spectralPSF_3D(rows,cols,:) = nan;
    s = spectralPSF_3D(:,:,1);  % get ready to use linear sampling.
    % check to see which ones are nan
    samp_ind_nonan = samp_ind(~isnan(s(samp_ind)));
    % use samp_ind_nonan for sampling from now on:
    spectralPSFrand = zeros(length(samp_ind_nonan),N3);
    
    % 4. use the indices without nan on all wavelength channels
    for m = 2:N3
        s = spectralPSF_3D(:,:,m);
        spectralPSFrand(:,m) = s(samp_ind_nonan); %use linear indexing for truly random
    end

end