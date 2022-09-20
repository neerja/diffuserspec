function [noisySignal] = generateNoise(signal,mu,stddev)
% This function's purpose is to generate a noise filter that simulates the
% noise present in optical imaging - specifically noise present in optical
% coherence tomography and spectroscopy

% SYNTAX:  [] = generateNoise(a,b,c)
%
% INPUT PARAMETERS:
%     - 
%
% OUTPUT PARAMETERS:
%     - 
%
% EXAMPLE: 
%     = generateNoise();
%
% Other m-files required: [none]
% Subfunctions: [none]
% MAT-files required: [none]
%
% See also: none

%--------------------------------------------------------------------------
%
%  Author:         Joe Malone - joseph.d.malone@vanderbilt.edu
%  Organization:    Bowden Biomedical Optics Lab, Vanderbilt University
%  Creation Date:   2021.02.15
%  Last Modified:   2021.02.15
%
%  Matlab Style Template v 1.0

%------------- BEGIN CODE --------------

noiseVec = mu + randn(size(signal)).*stddev;
noisySignal = signal+noiseVec;

end