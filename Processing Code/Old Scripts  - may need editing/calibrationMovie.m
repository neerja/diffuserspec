% calibrationMovie.m
% June 21st, 2022
% Neerja Aggarwal
% make movie of calibration frames
% assumes loaded spectralPSF_3D

figure;

writerObj = VideoWriter('out2','MPEG-4');
plotevery = 6
writerObj.FrameRate = 60/plotevery;
open(writerObj);
fid = figure;
imagesc(squeeze(spectralPSF_3D(:,:,1)));
%%
[m,n,l] = size(spectralPSF_3D);

for k = [1:plotevery:l]
    k
    pause(0.1);
    figure(fid);
    imagesc(squeeze(spectralPSF_3D(:,:,k)));
    axis off;
    name = [num2str(calibrationWavelengths(k)) ' nm'];
    title(name)
    frame =  getframe(gcf);
    writeVideo(writerObj,frame);
end
hold off

close(writerObj)

%% plot correlation matrix

% samp_line is the 2D matrix containing randomly smapled points per
% wavelength PSF.  samp_line size is (345x5476) (extra wavelength channel
% at end we can ignore)
% normalize each row;
[N1, N2] = size(samp_line);
for k = 1:N1
    v = samp_line(k,:);
    samp_norm(k,:) = v./norm(v);
end

corrmat = samp_norm*samp_norm';
figure
imagesc(corrmat)
xlabel('Wavelength Channel')
ylabel('Wavelength Channel')
title('Inner Product Matrix')
%%
figure(1);
plot(calibrationWavelengths_fit,corrmat(50,1:344));
hold on;
plot(calibrationWavelengths_fit,corrmat(300,1:344))
xlabel('Wavelength Channel')
ylabel('Inner Product')

thresh = 0.8
%zero out the top right section
min(abs(corrmat-thresh),1)
