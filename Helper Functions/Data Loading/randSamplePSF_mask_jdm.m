function [spectralPSFrand, samp_ind_nonan,mask] = randSamplePSF_mask_jdm(spectralPSF_3D,samplepercent,maskbox,pattern)

% samples an image randomly.  sets values inside center box of size
% maskbox (half width) to 0 and omits those elements in output.

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

    

    if pattern == 1 %middle omit

        mask = ones(N1,N2);
        % 2. get indices for mask
        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        mask(rows,cols) = 0;

    elseif pattern == 2 %corners

        mask = zeros(N1,N2);

        rows = [1:maskbox,N1-maskbox:N1];
        cols = [1:maskbox,N2-maskbox:N2];

        mask(rows,cols) = 1;

    elseif pattern == 3 % donut
        
        mask = zeros(N1,N2);

        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        rowsub = [floor(N1/2)-(maskbox/2):floor(N1/2)+(maskbox/2)];
        colsub = [floor(N2/2)-(maskbox/2):floor(N2/2)+(maskbox/2)];

        mask(rows,cols) = 1;
        mask(rowsub,colsub) = 0;

    elseif pattern == 4 %middle keep

        mask = zeros(N1,N2);
        % 2. get indices for mask
        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        mask(rows,cols) = 1;

    elseif pattern == 5 % invert donut
        
        mask = ones(N1,N2);

        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        rowsub = [floor(N1/2)-(maskbox/2):floor(N1/2)+(maskbox/2)];
        colsub = [floor(N2/2)-(maskbox/2):floor(N2/2)+(maskbox/2)];

        mask(rows,cols) = 0;
        mask(rowsub,colsub) = 1;

    elseif pattern == 6 % invert donut
        
        mask = ones(N1,N2);

        rows = [floor(N1/2)-maskbox:floor(N1/2)+maskbox];
        cols = [floor(N2/2)-maskbox:floor(N2/2)+maskbox];

        rowsub = [floor(N1/2)-(maskbox/4):floor(N1/2)+(maskbox/4)];
        colsub = [floor(N2/2)-(maskbox/4):floor(N2/2)+(maskbox/4)];

        mask(rows,cols) = 0;
        mask(rowsub,colsub) = 1;
    end

    % 3. set value to 0 inside A for mask
    spectralPSF_3D = spectralPSF_3D.*mask;

toc
    spectralPSF_3D = squeeze(reshape(spectralPSF_3D,N1*N2,1,N3));

    % check to see which ones are 0
    samp_ind_nonan = samp_ind(spectralPSF_3D(samp_ind) ~= 0);

    % take indexed nonNaN values for 2D PSF
    spectralPSFrand = spectralPSF_3D(samp_ind_nonan,:);
toc
end