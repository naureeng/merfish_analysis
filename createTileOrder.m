function [tileOrder] = createTileOrder(tileSize,option)

% function that takes in the size n x m of the desired tiles and creates a
% corresponding matrix according to the option specified

nRows  = tileSize(1);
nCols  = tileSize(2);
nTiles = nRows*nCols;

switch option
    case 'snake'
        to_flat     = 1:nTiles;
        tileOrder   = reshape(to_flat,tileSize);    % this creates a top down tile order matrix
                                                % to get it into snake format, we need to flip the first and then every second column up down
        nFlipCols   = round(nCols/2);           % this tells us how many columns need to be flipped up down
        % figure out which is the last column we need to flip
        if rem(nCols,2) == 0                    % if this is true, the number of columns is even, which means the second to last column is the last column to be flipped                    
            lastFlipCol = nCols - 1;
        else                                    % otherwise the last column to be flipped is the last column
            lastFlipCol = nCols;
        end
        flipColumns = linspace(1,lastFlipCol,nFlipCols);
        % flip the columns
        tileOrder(:,flipColumns)   = flipud(tileOrder(:,flipColumns));
        
    case 'linear'
        tileOrder = 1:nTiles;
        
    otherwise
        error('unknown tile order option');
end

end