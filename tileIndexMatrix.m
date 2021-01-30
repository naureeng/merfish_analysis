function [tileIndices] = tileIndexMatrix(tileSize)

% function that calculates the indices within the big tiled image for the
% individual tiles to go into
% tileSize is the m x n matrix for how many tiles go into one row and one column, respectively, in the big image

% tileIndices is a m x n cell that within each cell, contains the row and
% column indices that the corresponding tile will occupy

rowPix = 2048;
colPix = 2048;
pixelDimensions = [rowPix colPix];

nRows = tileSize(1);
nCols = tileSize(2);
nTiles = nRows*nCols;

tileMatrix = reshape(1:nTiles,tileSize);

nPixRow = nRows*rowPix;
nPixCol = nCols*colPix;

tileIndices = {};

for tc = 1:nCols
    for tr = 1:nRows
        
        tileIndices{tr,tc}(1,:) = [1:rowPix]+((tr-1)*rowPix);
        tileIndices{tr,tc}(2,:) = [1:colPix]+((tc-1)*colPix);
        
    end
    
end

end