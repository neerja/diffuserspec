function [diffuserCoefficients] = MakeDiffuserMatrix(nPixels)
% MakeDiffuserMatrix - Create a matrix of random coefficients for wavelength
% mapping to mimic action of a diffuser on OCT spectral data.
% Assumes no loss (all weights sum to 100% of input light for each
% wavelength). That is, coefficients add to 1 for each row. Assumes no
% saturation - that is, columns do not necessarily add up to 1.
%
% SYNTAX:  [dc] = MakeDiffuserMatrix(nPixels)
%
% INPUT PARAMETERS:
%    nPixels - number of pixels in the spectrometer
%
% OUTPUT PARAMETERS:
%    diffuserCoefficients - square matrix with coefficients for each
%       wavelength (row = pixel#, column = wavelength#)
%
% EXAMPLE: 
%    dc = MakeDiffuserMatrix(1024);
%
% Other m-files required: [none]
% Subfunctions: [none]
% MAT-files required: [none]
%
% See also: none

%--------------------------------------------------------------------------
%
%  Author:          Audrey Bowden - a.bowden@vanderbilt.edu
%  Organization:    Bowden Biomedical Optics Lab, Vanderbilt University
%  Creation Date:   2020.12.16
%  Last Modified:   2020.12.16
%
%  Matlab Style Template v 1.0


%------------- BEGIN CODE --------------

%setup diffuser matrix
diffuserCoefficients = rand(nPixels);
diffuserCoefficients = diffuserCoefficients./...
    repmat(sum(diffuserCoefficients,2),1,nPixels);

%------------- END OF CODE --------------
end