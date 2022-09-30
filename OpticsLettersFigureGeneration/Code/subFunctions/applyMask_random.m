function [SSTM_2D,samp_ind] = applyMask_random(SSTM,maskCoordinates,samplepercent)

    numCoords = size(maskCoordinates,1);
    numsamples = floor(numCoords*(samplepercent/100));
    samp_ind = randi(numCoords,[numsamples,1]);

    SSTM_2D = zeros(numsamples,size(SSTM,3));

    for coordinate = 1:numsamples
        SSTM_2D(coordinate,:) = SSTM(maskCoordinates(samp_ind(coordinate,1),2)...
            ,maskCoordinates(samp_ind(coordinate,1),1),:);

        if mod(coordinate,5000) == 0
            disp([num2str(coordinate)])
        else
        end

    end
end