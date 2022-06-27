%load transfermatrix
transferMatrix_3D = tifRead('220grit_transferMatrix_3D_normalized_8bit.tif');

%load multipeak data
multiPeak = tifRead('multiPeak_diffuseImage.tif');

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

%extract line from diffuse data for recon
b = multiPeak(lineNum,:)';

%apply SVD threshold
thresh = .0125;

%perform truncated SVD recon
[recon,~,Dprime,Dprime_thresh] = TruncSVD(A_crop,b,thresh);

wave =linspace(calibrationWavelengths(1),calibrationWavelengths(end),173);

figure;
plot(wave, abs(recon))
xlabel('wavelength (nm)')


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

function [x,Ainv,Dprime,Dprime_thresh] = TruncSVD(A,b,thresh)

    % A = calibrated diffuser transform matrix
    [U,D,V] = svd(A);

    % D = diagnol matrix containing singular values
    Dprime = (1./D)';
    Dprime(Dprime == inf) = 0;
    Dprime_thresh = Dprime;
    Dprime_thresh(Dprime_thresh > thresh) = 0;

    % Solve for inverse of diffuser transform matrix
    Ainv = V*Dprime_thresh*U';

    %perform inversion, x = reconstructed spectrum
    x = Ainv*b;
    
end