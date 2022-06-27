%% process_220grit
% June 7th, 2022 Data
% Neerja Aggarwal
% Purpose: Remake the calibration matrix.  


%% import file

cd JoeData_07062022/

%import the 3D calibration file using special tif Read function
%(implemented below)
transferMatrix_3D = tifRead('220grit_transferMatrix_uint16.tif');
load calibrationWavelengths_polyfit.mat

%% check power of 3D matrix
tm3Dsub = transferMatrix_3D-bkgd_50ms;
%%
pwr = squeeze(sum(tm3Dsub,[1,2]));

figure
plot(wvcal,pwr(1:344))
xlabel('Wavelength')
ylabel('Intensity')
%%
pwrnorm = pwr(1:344)./wvint;
figure
plot(pwrnorm)

%% normalize the 2dsub using pwr from 3d. 
pwr = pwr./max(pwr);

tm2Dsub = transferMatrix_2Dsub./pwr';

pwr2D = sum(tm2Dsub,1);

figure
plot(wvcal,pwr2D(1:344))
ylim([0, 5.5e6])
xlabel('Wavelength')
ylabel('Intensity')
%% make 2D matrix
% select line number to use for 2D transfer matrix
lineNum = 900;
transferMatrix_2D = squeeze(transferMatrix_3D(lineNum,:,:));

% check for saturation
maxval = max(max(transferMatrix_2D))
% we're good

%% subtract background and normalize using the spf2 files

bkgd_50ms = mean(tifRead('220grit_transferMatrix_background_50ms.tif'),3);
bkgd_50ms_2D = bkgd_50ms(lineNum,:)';

transferMatrix_2Dsub = transferMatrix_2D-bkgd_50ms_2D;

pwr = sum(transferMatrix_2Dsub,1);
plot(calibrationWavelengths_polyfit,pwr(1:344))
xlabel('Wavelength (nm)')
ylabel('Intensity')

%%

% createSourceSpectrum.m
% Neerja Aggarwal
% June 8th, 2022
% Purpose: import the spf2 spectrum. find peak (using either max function
% or findpeak function).  store the wavelength and intensity. 

cd transferMatrix_gtWavelengths

files = dir('*.spf2')

N = length(files)

wvcal = zeros(N,1)
wvint = zeros(N,1)

for k1 = 1:N
    filename = files(k1).name;
    % import
    spectrum = readSPF2(filename);
    % get max val
    [val,idx] = max(spectrum(:,2));
    wvcal(k1) = spectrum(idx,1);
    wvint(k1) = val;
end
cd ..
wvint = wvint./max(wvint)

figure
plot(wvcal,wvint)
xlabel('Wavelength (nm)')
ylabel('Intensity')

transferMatrix_norm = transferMatrix_2Dsub(:,1:344)./(wvint');
%
figure
pwr = sum(transferMatrix_norm,1);
plot(calibrationWavelengths_polyfit,pwr)
xlabel('Wavelength (nm)')
ylabel('Intensity')
ylim([0,12e6])

%% (SKIP) normalize by the incident power using broadband_gt.spf2
% load the ground (incident spectrum)
broadband_gt = readSPF2('broadband_gt.spf2')
bbgt_int = interp1(broadband_gt(:,1),broadband_gt(:,2),calibrationWavelengths_polyfit);

figure
plot(calibrationWavelengths_polyfit,bbgt_int)
% use this graph to brush, flatten 850 nm peak, and make sourceSpectrum
% array normalize the source spectrum
sourceSpectrum(:,2) = sourceSpectrum(:,2)./max(sourceSpectrum(:,2));
% normalize the transfermatrix
transferMatrix_norm = transferMatrix_2D(:,1:344)./(sourceSpectrum(:,2)');

% now we have our calibration file.  Matlab automatically converted the
% uint1t6->double earlier. 
%%
% Now we need to bin our calibration file so that we have ~2 nm bins. This
% will give us 42 wavelength channels which is reasonable: 870-784 = 86
% nm/2 = 43 -1
% channels
wavevec = [785:2:869];

% find the element in calbirationWavelengths_polyfit that matches the
% current wv value:
wvst = wavevec(1)
[~, wvind] = min(abs(calibrationWavelengths_polyfit-wvst))
calibrationMatrix = zeros(2560,42);
for k1 = 1:length(wavevec)-1
    wvend = wavevec(k1+1);
    [~, wvind2] = min(abs(calibrationWavelengths_polyfit-wvend));
    wvind2
    
    % go from wvind to wvind2-1
    cut = transferMatrix_norm(:,wvind:wvind2);
    cutsum = sum(cut,2);
    calibrationMatrix(:,k1) = cutsum;
    
    wvst = wvend;
    wvind = wvind2;
end


%% try another way - sum every 8 frames
sumFrames = 8
[N1,N2] = size(transferMatrix_norm)

calibrationMatrix = zeros(N1,floor(N2/sumFrames));
wvvec2 = zeros(floor(N2/sumFrames),1);

for k1 = 1:length(wvvec2)
wvvec2(k1) = calibrationWavelengths_polyfit((k1-1)*sumFrames+1);
cutsum = mean(transferMatrix_norm(:,(k1-1)*sumFrames+1:k1*sumFrames),2);
calibrationMatrix(:,k1) = cutsum;
end

figure
imagesc(calibrationMatrix)

figure
plot(wvvec2,'*-')
ylabel('Wavelength (nm)')
xlabel('Index')

save calibrationMatrix1.mat wvvec2 calibrationMatrix lineNum

%% measurements using CalibrationMatrix2:

uniformsource = ones(length(wvvec2),1);
uniform_meas = calibrationMatrix2*uniformsource;
figure
plot(uniform_meas)


%%
%normalize each vector first:
tmnorm = zeros(N1, N2);
for k1 = 1:N
    tmnorm(:,k1) = transferMatrix_norm(:,k1)./norm(transferMatrix_norm(:,k1));
end

coherence = tmnorm'*tmnorm;
figure
imagesc(coherence)
colorbar()

%%
idx1 = 30;
idx2 = N-idx1;
figure
plot(calibrationWavelengths_polyfit,coherence(idx1,:))
hold on
plot(calibrationWavelengths_polyfit,coherence(idx2,:))
xlabel('Wavelength (nm)')
ylabel('Intensity')

%%
function [stack] = tifRead(filename,maxFrame,x,y)
%
%
% filename = name of .tiff or multitiff you want to read
% x = [start,end] - pixel height boundry; i.e. [400,800] will crop and read in only 
%       pixel range 400:800
% y = [start,end] - pixel width boundry
%
    
    info = imfinfo(filename);
    nframes = size(info,1);
    if exist('maxFrame','var')
    else 
        maxFrame = nframes;
    end
    
    if exist('x','var') && exist('y','var')
    
    elseif exist('x','var') == 0 && exist('y','var')
        x = [1,info(1).Height];
        
    elseif exist('x','var') && exist('y','var') == 0 
        y = [1,info(1).Width];
        
    else
        x = [1,info(1).Height];
        y = [1,info(1).Width];
    end
    
    im = imread(filename,1,'PixelRegion',{x,y});
    stack = zeros(size(im,1),size(im,2),maxFrame);
    stack(:,:,1) = im;
        
    for n = 2:maxFrame
    
        stack(:,:,n) = imread(filename,n,'PixelRegion',{x,y});
    
    end

end

