%load transfermatrix
transferMatrix_3D = tifRead('220grit_transferMatrix_3D_normalized_8bit.tif');

%load multipeak data

multiPeak = tifRead('multiPeak_diffuseImage.tif');
broadband = tifRead('broadband_diffuse_#0001.tif');

broadband_avg = mean(broadband,3);

%load calibrated ground truth wavelength data

calibrationWavelengths = load('calibrationWavelengths.mat');
calibrationWavelengths = calibrationWavelengths.calibrationWavelengths;

%% Process and reconstruct multipeak data

% subsample transfer matrix by a factor of 2;
transferMatrix_3D_crop = transferMatrix_3D(:,:,1:2:end);

% select line number to use for 2D transfer matrix
lineNum = 900;

%extract line data from transfer matrix
A_crop = squeeze(transferMatrix_3D_crop(lineNum,:,:));
Afull_crop = squeeze(transferMatrix_3D(lineNum,:,:));

%extract line from diffuse data for recon
b = broadband_avg(lineNum,:)';

%apply SVD threshold
thresh = .0125; % 
% thresh = .0001; % 
% thresh = .001; % 

%perform truncated SVD recon
%[recon,Alr,Ainv,Dprime,Dprime_thresh] = TruncSVD(A_crop,b,thresh);
[recon,Alr,Ainv,Dprime,Dprime_thresh] = TruncSVD(Afull_crop,b,thresh); % use the full A

wave =linspace(calibrationWavelengths(1),calibrationWavelengths(end),173);

figure;
%plot(wave, abs(recon))
plot(calibrationWavelengths,abs(recon(1:end-1)))
xlabel('wavelength (nm)')

%% try to see the ground truth

broadband_gt = readSPF2('broadband_gt.spf2')

figure
plot(broadband_gt(:,1),broadband_gt(:,2))

% try to pass it through the forward model and see what the measurement
% should look like.

% interpolate the spectrum
bbgt_int = interp1(broadband_gt(:,1),broadband_gt(:,2),calibrationWavelengths_polyfit);

% plot
figure
plot(calibrationWavelengths,bbgt_int)
xlabel('Wavelength (nm)')
ylabel('Intensity')

%forward model (unnormalized)
bsim = Afull_crop(:,1:344)*bbgt_int;

%% try to solve from simulated measurement after adding a little bit of noise
uniform = sum(Afull_crop,2)
cdif  = uniform/max(uniform)-bsim/max(bsim)
figure
plot(cdif)


%% look at Afull_crop matrix

figure
imagesc(Afull_crop)
xlabel('Wavelength Channel')
ylabel('Line 900')

%% subfunctions

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

function [x,Alr,Ainv,Dprime,Dprime_thresh] = TruncSVD(A,b,thresh)

    % A = calibrated diffuser transform matrix
    [U,D,V] = svd(A);

    % D = diagnol matrix containing singular values
    Dprime = (1./D)';
    Dprime(Dprime == inf) = 0;
    Dprime_thresh = Dprime;
    Dprime_thresh(Dprime_thresh > thresh) = 0;
    Dthresh = D;
    Dthresh(Dthresh< (1./thresh)) = 0;

    % Solve for inverse of diffuser transform matrix
    Ainv = V*Dprime_thresh*U';
    Alr = U*Dthresh*V';

    %perform inversion, x = reconstructed spectrum
    x = Ainv*b;
    
end