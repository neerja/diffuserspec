
function [x,Alr,Ainv,Dprime,Dprime_thresh] = TruncSVD(A,b,thresh)

    % A = calibrated diffuser transform matrix
    [U,D,V] = svd(A);

    % D = diagnol matrix containing singular values
    Dprime = (1./D)';
    Dprime(Dprime == inf) = 0;
    Dprime_thresh = Dprime;
    
    % set any singular values (inverted) above threshold to zero
    Dprime_thresh(Dprime_thresh > thresh) = 0;

    % Solve for inverse of diffuser transform matrix
    Ainv = V*Dprime_thresh*U';
    
    % use this code to also output the low rank A:
    Dthresh = D;
    % set singular values below 1/threshold to zero (similar effect)
    Dthresh(Dthresh< (1./thresh)) = 0;
    Alr = U*Dthresh*V';

    %perform inversion, x = reconstructed spectrum
    x = Ainv*b;
    
end