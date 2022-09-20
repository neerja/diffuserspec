function [mask,maskCoordinates] = makeCircularMask(data,x,y,radius)

xyrVec = [x,y,radius];

[xSize,ySize,~] = size(data);

[X,Y] = meshgrid(1:xSize, 1:ySize);  %create pixel coordinates

mask = hypot(X - xyrVec(1), Y - xyrVec(2));
mask = any(mask<= xyrVec(3), 3);

[xLoc,yLoc] = find(mask);
maskCoordinates = cat(2,xLoc,yLoc);

end