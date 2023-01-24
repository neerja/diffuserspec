% Purpose: main reconstruction algorithm.  uses truncated SVD to reconstruct spectra
% from measurement of speckle.  but uses a gaussian filter instead of box 
% to filter the singular values. 
%Inputs: A - 2D SSTM, obtained from calibration
%        b - input diffuse spectrum to be reconstructed
%        sigma - gaussian filter width (affects rolloff)

%Outputs: 
%       recon = reconstructed input spectrum (main output)
%       dvec_inv = inverse singular values (for plotting purposes)
%       dvec_invf = filtered inverse singular values (for plotting purposes)
%       gaussian_filter = filter for inverse singular values (for plotting purposes)

function [recon,dvec_inv,dvec_invf,Ainvf] = gaussSVD(A,b,sigma)


% initialize size of input matrix A
    [numPixels,numWavelengths] = size(A);

% perform SVD on SSTM
%     [U,D,V] = svd(A);
    [U,D,V] = svd(A,'econ');

% create singular value vector by obtaining diagnol values of D
    dvec = diag(D);
    indexvec = 1:length(dvec);

mu = 1; % (center for gauss filter)
a = 1; % scaling coefficient

% calculate normalized gaussian filter
    gaussian_filter = exp(-(indexvec-mu).^2/(2*sigma^2)) / (sigma*sqrt(2*pi));
    gaussian_filter = a*gaussian_filter/max(gaussian_filter); %normalize and scale

% invert singular value vector and multiply by gaussian filter
    dvec_inv = 1./dvec;
    dvec_invf = dvec_inv.*gaussian_filter';

% make diagnol matrix from inverse of D
    Df_inv = diag(dvec_invf);
%     Z = zeros(numWavelengths,abs(numPixels-numWavelengths));
%     Df_inv = [Df_inv Z]; 

% solve for inverse of A matrix
    Ainvf =V*Df_inv*U';
    
% perform reconstruction  
    recon = Ainvf*b;


end