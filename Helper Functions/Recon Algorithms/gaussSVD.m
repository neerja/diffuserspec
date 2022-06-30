function x = gaussSVD(A,b,sigma)

[m,n] = size(A);
[U,D,V] = svd(A);
dvec = diag(D);

if isempty(sigma)
    sigma = 30 ; % (std dev for gauss filter)
end

indexvec = 1:length(dvec);
mu = 1; % (center for gauss filter)
a = 1; % scaling coefficient
gaussian_filter = exp(-(indexvec-mu).^2/(2*sigma^2)) / (sigma*sqrt(2*pi));
gaussian_filter = a*gaussian_filter/max(gaussian_filter); %normalize and scale

dvec_inv = 1./dvec;
dvec_invf = dvec_inv.*gaussian_filter';

% make inverse and do recon
Df_inv = diag(dvec_invf);
Z = zeros(n,abs(m-n));
Df_inv = [Df_inv Z]; 

Ainvf =V*Df_inv*U';

% b = spectrumForRecon_sampled; 
x = Ainvf*b;
x = x./max(x);

end