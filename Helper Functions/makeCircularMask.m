% Purpose: create a matrix of mask coordinates 
% Inputs:
% data = this is the SSTM (3D matrix)  
% innerRad = scalar
% Set innerRad = 0 for full circular mask
% Outputs:
% mask = 2D image containing 1's and 0's at specified locations according to the mask
% maskCoordinates = nx2 array containing pixel coordinates for 1's. (used in applyMask_random.m)

function [mask,maskCoordinates] = makeCircularMask(data,innerRad,outerRad)

[xSize,ySize,~] = size(data);

xyrInnerVec = [xSize/2,ySize/2,innerRad];
xyrOuterVec = [xSize/2,ySize/2,outerRad];

[X,Y] = meshgrid(1:xSize, 1:ySize);  %create pixel coordinates

maskInner = hypot(X - xyrInnerVec(1), Y - xyrInnerVec(2));
maskInner = any(maskInner<= xyrInnerVec(3), 3);

if outerRad == 0
    
    maskOuter = true(ySize,xSize);
    
elseif outerRad
    
    maskOuter = hypot(X - xyrOuterVec(1), Y - xyrOuterVec(2));
    maskOuter = any(maskOuter<= xyrOuterVec(3), 3);
    
end

mask = maskOuter-maskInner;

[xLoc,yLoc] = find(mask);
maskCoordinates = cat(2,xLoc,yLoc);

end