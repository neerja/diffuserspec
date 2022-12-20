function [mask,maskCoordinates] = makeCircularMask(data,innerRad,outerRad)

% Set innerRad = 0 for full circular mask

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