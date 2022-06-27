%% process_broadbanddiffuse.m
% Neerja Aggarwal
% June 7th, 2022

% load broadband diffuse measurement file.  Make HDR measurement. load
% calibration file.  Use line number data to do pseudoinverse. 

load calibrationMatrix2.mat
bbdiff1 = mean(tifRead('broadband_diffuse_#0001_uint16.tif'),3); 
bbdiff1_line = bbdiff1(lineNum,:);
figure
plot(bbdiff1_line)

bbdiff2 = mean(tifRead('broadband_diffuse_#0002_uint16.tif'),3); 
bbdiff2_line = bbdiff2(lineNum,:);
figure
plot(bbdiff2_line)

bbdiff3 = mean(tifRead('broadband_diffuse_#0003_uint16.tif'),3); 
bbdiff3_line = bbdiff3(lineNum,:);
figure
plot(bbdiff3_line)
% FRAME 3 looks good lets use this.

%% pseudoinv the measurement: 
A = calibrationMatrix2;
b = bbdiff3_line';

%least squares solution:
xest = A\b

% residual


figure
plot(wvvec2,xest)
xlabel('Wavelength (nm)')
ylabel('Intensity')


%% look at the correlation matrix:

cmatrix = A'*A;

figure
imagesc(cmatrix)

%% look at the raw calibration speckles:

figure
plot(A(:,1:5:end))

pwr = mean(A,1)
figure
plot(wvvec2,pwr)
xlabel('Wavelength (nm)')
ylabel('Intensity')

%% look at the original transfer matrix:

pwr = mean(transferMatrix_2D,1)
figure
plot(calibrationWavelengths_polyfit,pwr(1:344))
xlabel('Wavelength (nm)')
ylabel('Intensity')

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