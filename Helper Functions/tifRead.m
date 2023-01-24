% Purpose: read tif file containing calibration matrix data.
% Inputs: 
% filename = name of .tiff or multitiff you want to read
% x = [start,end] - pixel height boundary; i.e. [400,800] will crop and read in only 
%       pixel range 400:800
% y = [start,end] - pixel width boundary; i.e. [400,800] will crop and read in only 
%       pixel range 400:800
% Outputs:
% stack = 3D array containing calibration matrix with axes: x,y,lambda

function [stack] = tifRead(filename,maxFrame,x,y)

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