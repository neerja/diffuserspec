% loadTestMeasurement.m
% Function inputs: filename of the measurement .tif
% filename of the background tif, sampxy is the sampling indices
% sampxy and sampfac are optional to sample specific pixels. 
function [spectrumForRecon_sampled] = loadTestMeasurement(filename,bgfilename,samp_xy,sampfac)

spectrumForRecon = tifRead(filename);
bgSpectrumForRecon = tifRead(bgfilename);
% average along 3rd dim
spectrumForRecon = mean(spectrumForRecon,3);
bgSpectrumForRecon = mean(bgSpectrumForRecon,3);
spectrumForRecon = abs(spectrumForRecon - bgSpectrumForRecon);
if isempty(samp_xy)
    % if no sampling, then return entire frame
    spectrumForRecon_sampled = spectrumForRecon;
else
    spectrumForRecon_sampled = spectrumForRecon(samp_xy(:,1),samp_xy(:,2));
    spectrumForRecon_sampled = reshape(spectrumForRecon_sampled,sampfac^2,1);
end
end