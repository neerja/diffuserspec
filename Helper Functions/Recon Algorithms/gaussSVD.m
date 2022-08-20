function x = gaussSVD(A,b,sigma)

%Inputs: A - 2D SSTM, obtained from calibration
%        b - input diffuse spectrum to be reconstructed
%        sigma - gaussian filter kernal for thresholding

%Outputs: x - reconstructed input spectrum

% initialize size of input matrix A
[m,n] = size(A);

% perform SVD on SSTM
[U,D,V] = svd(A);

% create singular value vector by obtaining diagnol values of D
dvec = diag(D);

% if no sigma is assigned, sigma default = 30
if isempty(sigma)
    sigma = 30 ; % (std dev for gauss filter)
end

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
Z = zeros(n,abs(m-n));
Df_inv = [Df_inv Z]; 

% solve for inverse of A matrix
Ainvf =V*Df_inv*U';

% perform reconstruction  
x = Ainvf*b;

%x = x./max(x);

end