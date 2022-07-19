function [spectralPSF_2D, samp_xy,sampfac] = randSamplePSF_mask(spectralPSF_3D,samplepercent,maskSize)

%This function doesn't do the right things

[N1,N2,N3] = size(spectralPSF_3D);
sampfac = round(sqrt(N1*N2*(samplepercent/100)));

N1_mask = N1-maskSize;
N2_mask = N2-maskSize;

samp_xy = zeros(sampfac,2);
samp_xy(:,1) = ceil(N1_mask*rand(sampfac,1))+maskSize;
samp_xy(:,2) = ceil(N2_mask*rand(sampfac,1))+maskSize;

samp_full = spectralPSF_3D(samp_xy(:,1),samp_xy(:,2),:);
spectralPSF_2D = reshape(samp_full,[sampfac^2,N3]);

end