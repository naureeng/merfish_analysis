function createTiledImage(expName,resRound,whichGenes,stackNumbers,tileOrder)

% function that tiles stacks together & overlays the gene localisations
% onto the corresponding DAPI images

% inputs
% expName       name of the experiment we want to look at
% resRound      name of the round of storm analysis results we want to look at
% whichGenes    a list of genes that we want to visualise
%               must be in the same format as the MERFISH analysis results
% stackNumbers  which tiles we want to put together

% optional
% tileOrder     a matrix that tells us in which order the tiles need to go

%%
% add an option to run the createTileOrder function here if it doesn't exist yet
if nargin<5
    disp('tile order not provided.');
    tileSize        = input('provide tilesize (m rows n columns): ');       % this needs to be 2 numbers inside square brackets []
    tilingOption    = input('provide tiling option (eg snake): ');          % this needs to be given as a string, so like this: 'snake'
    [tileOrder] = createTileOrder(tileSize,tilingOption);
end

% check that inputs make sense
if length(stackNumbers) ~= size(tileOrder,1)*size(tileOrder,2)
    error('stack numbers and tile order do not match up');
end

nTilesRow = size(tileOrder,1);          % how many tiles will be in a row
nTilesCol = size(tileOrder,2);          % how many tiles will be in a column
nTiles    = nTilesRow*nTilesCol;        % total number of tiles

% set up the paths where the data lives
geneLocsPath    = fullfile('D:\MERFISH results',expName,resRound,'geneLocs');          % path for the folder where the XY coordinates for the genes live
dapiFilesPath   = fullfile('D:\processed',expName);

dapiFileName    = 'bit1_dapi.tif';      % for now, use the first round dapi image as the background image

pixelDimensions = [2048 2048];          % size of the original tif files in pixels
nPix            = pixelDimensions(1)*pixelDimensions(2);    % total number if pixels in a tile
nGenes          = length(whichGenes);   % number of genes provided
% determine final matrix size
finalImage = zeros(pixelDimensions(1)*nTilesRow,pixelDimensions(2)*nTilesCol);

% calling the function that calculates the indices within the final image that the individual tiles will occupy
inds    = tileIndexMatrix([nTilesRow nTilesCol]);
inds1D  = reshape(inds,nTiles,1);                   % reshape it so it's more like a vector, which makes things easier in the later loops


%% create a tiled dapi image first

% check if the tiled dapi image exists already, else create it
if ~exist(fullfile('D:\MERFISH results',expName,resRound,...
        strcat('tiledDapiImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.fig')));
    
    fullDapiImage  = finalImage;
    
    for tt = 1:nTiles
        
        tileNum     = tileOrder(tt);            % this tells us the how manieth tile we are putting in the current position
        stackNum    = stackNumbers(tileNum);    % this tells us what stack is supposed to go into this position
        thisDapiFile = fullfile(dapiFilesPath,...
            strcat('zscan_',num2str(stackNum.','%02d')),dapiFileName);  % determine what the filename of the stack we want to load is
        dapiIm = imread(thisDapiFile);          % load the image
        
        whichInds = inds1D{tt,1};               % this tells us the indices within the full image that this stack is going to occupy
        
        fullDapiImage(whichInds(1,:),whichInds(2,:)) = im2double(dapiIm);  % add the image to its position in the big tiled image
        
    end
    
    % plot the tiled dapi image
    fd = figure;
    imagesc(fullDapiImage), axis off; axis equal
    title([expName ' stacks ' num2str(stackNumbers(1)) ' to ' num2str(stackNumbers(end))]);
    
    %save tif file in merfish results folder
    imwrite(fullDapiImage,fullfile('D:\MERFISH results',expName,resRound,...
        strcat('tiledDapiImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.tif')));
    
    savefig(fd,fullfile('D:\MERFISH results',expName,resRound,...
        strcat('tiledDapiImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.fig')));
    
else
    
    fullDapiImage = imread(fullfile('D:\MERFISH results',expName,resRound,...
        strcat('tiledDapiImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.tif')));
    
    fullDapiImage = im2double(fullDapiImage);
    
end

%% load gene files

% go through the genes we want to look at in a loop
for g = 1:nGenes
    thisGene    = whichGenes{g};
    
    if ~exist(fullfile(geneLocsPath,thisGene,...
            strcat(thisGene,'_tiledImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.fig')));
        
        geneImage   = finalImage;
        
        for tt = 1:nTiles
            
            tileNum     = tileOrder(tt);            % this tells us the how manieth tile we are putting in the current position
            stackNum    = stackNumbers(tileNum);    % this tells us what stack is supposed to go into this position
            thisGeneFile = fullfile(geneLocsPath,thisGene,...
                strcat('zscan_',num2str(stackNum.','%02d'),'_',thisGene,'_XYcoordinates.mat'));
            if ~exist(thisGeneFile)                 % for the very first experiment, the filename is slightly different
                thisGeneFile = fullfile(geneLocsPath,thisGene,...
                    strcat('zstack_zscan_',num2str(stackNum.','%02d'),'_',thisGene,'_XYcoordinates.mat'));
            end
            
            % in some stacks, there may be no expression of the gene of interest, which is why we need to add this if statement, otherwise there will be an error for trying to load an inexistant file
            if exist(thisGeneFile)                  load(thisGeneFile);
                
                whichInds = inds1D{tt,1};           % this tells us the indices within the full image that this stack is going to occupy
                
                geneImage(whichInds(1,:),whichInds(2,:)) = imrotate(flipud(stackIm),-90);      % add the image to its position in the big tiled image
            end
            
        end
        
        % add color to the pixels around the main spot
        paddedGeneIm = zeros(pixelDimensions(1)*nTilesRow+2,pixelDimensions(2)*nTilesCol+2,9);
        paddedGeneIm(2:end-1,2:end-1,1) = geneImage;        % puts the original gene image into the centre
        
        paddedGeneIm(1:end-2,1:end-2,2) = geneImage*0.8;    % fills in the pixels on the diagonal left upper corner with 0.8
        paddedGeneIm(1:end-2,2:end-1,3) = geneImage*0.8;    % fills in the pixel on top with 0.8
        paddedGeneIm(1:end-2,3:end,4)   = geneImage*0.8;    % fills in the pixel on the diagonal right upper corner with 0.8
        
        paddedGeneIm(2:end-1,1:end-2,5) = geneImage*0.8;    % fills in the pixel to the left with 0.8
        paddedGeneIm(2:end-1,3:end,6)   = geneImage*0.8;    % fills in the pixel to the right with 0.8
        
        paddedGeneIm(3:end,1:end-2,7)   = geneImage*0.8;    % fills in the pixel on the diagonal left bottom corner with 0.8
        paddedGeneIm(3:end,2:end-1,8)   = geneImage*0.8;    % fills in the pixel below with 0.8
        paddedGeneIm(3:end,3:end,9)     = geneImage*0.8;    % fills in the pixel on the diagonal right bottom corner with 0.8
        
        paddedGeneIm = sum(paddedGeneIm,3);
        
        % uncomment to show this
        %     % plot spatial gene distribution by itself, with only the 'true' spots coloured
        %     figure;
        %     imagesc(geneImage), axis off; axis equal; colormap(gray); caxis([0 1]);
        %     title({[expName ' stacks ' num2str(stackNumbers(1)) ' to ' num2str(stackNumbers(end))];...
        %         [thisGene ' expression']});
        
        % plot the gene distribution with the neighbouring pixels shaded
        fg = figure;
        imagesc(paddedGeneIm), axis off; axis equal; colormap(gray); caxis([0 1]);
        title({[expName ' stacks ' num2str(stackNumbers(1)) ' to ' num2str(stackNumbers(end))];...
            [thisGene ' expression']});
        % save figure as tif
        imwrite(paddedGeneIm,fullfile(geneLocsPath,thisGene,...
            strcat(thisGene,'_tiledImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.tif')));
        savefig(fg,fullfile(geneLocsPath,thisGene,...
            strcat(thisGene,'_tiledImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.fig')));
        
        close all;
        
    else
        paddedGeneIm = imread(fullfile(geneLocsPath,thisGene,...
            strcat(thisGene,'_tiledImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.tif')));
        
    end
    
    if ~exist(fullfile(geneLocsPath,thisGene,...
            strcat(thisGene,'_tiledImage_withDAPI_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.fig')));
        
        % lastly overlay with dapi
        paddedDapi = zeros(size(paddedGeneIm))*NaN;
        paddedDapi(2:end-1,2:end-1) = fullDapiImage;
        figure;
        imshow(paddedDapi)      % show the dapi in black and white in the background
        caxis([0 0.4]);         % change these values to change how bright dapi looks
        hold on                 % this keeps the handle of the figure available for adding to it
        % make truecolor all green image for the gene spots to be overlaid on
        green = cat(3,zeros(size(paddedDapi)),ones(size(paddedDapi)),zeros(size(paddedDapi)));
        h = imshow(green);
        hold off
        set(h,'AlphaData',paddedGeneIm);
        title({[expName ' stacks ' num2str(stackNumbers(1)) ' to ' num2str(stackNumbers(end))];...
            ['round 1 DAPI with ' thisGene ' expression']});
        
        saveas(gcf,fullfile(geneLocsPath,thisGene,...
            strcat(thisGene,'_tiledImage_withDAPI_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)))),'pdf');
        
%         print(fullfile(geneLocsPath,thisGene,...
%             strcat(thisGene,'_tiledImage_withDAPI_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)))),'-dtiffn');
        
    end
    
    close all;
    
end

end

