function [spectralPSFrand, samp_ind_nonan] = randSamplePSF_mask_jdm(spectralPSF_3D,samplepercent,maskbox)

% samples an image randomly.  sets values inside center box of size
% maskbox (half width) to NaN and omits those elements in output.
tic
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
toc
    spectralPSF_3D = squeeze(reshape(spectralPSF_3D,[],1,N3));

    % check to see which ones are nan
    samp_ind_nonan = samp_ind(~isnan(spectralPSF_3D(samp_ind)));

    % take indexed nonNaN values for 2D PSF
    spectralPSFrand = spectralPSF_3D(samp_ind_nonan,:);
toc
end