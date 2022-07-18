function [spectralPSF_2D, samp_xy,sampfac] = randSamplePSF(spectralPSF_3D,samplepercent)


[N1,N2,N3] = size(spectralPSF_3D);
sampfac = round(sqrt(N1*N2*(samplepercent/100)));

samp_xy = zeros(sampfac,2);
samp_xy(:,1) = ceil(N1*rand(sampfac,1));
samp_xy(:,2) = ceil(N2*rand(sampfac,1));

samp_full = spectralPSF_3D(samp_xy(:,1),samp_xy(:,2),:);
spectralPSF_2D = reshape(samp_full,[sampfac^2,N3]);

end