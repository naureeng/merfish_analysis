function multiGeneOverlays(expName,resRound,whichGenes,stackNumbers,overlayColours)

% function that overlays several gene expression results over DAPI

% the createTiledImage function MUST have been run before this so that the
% corresponding tiled DAPI & individual tiled gene expression images exist

% inputs
% expName       name of the experiment we want to look at
% resRound      name of the round of storm analysis results we want to look at
% whichGenes    a list of genes that we want to visualise
%               must be in the same format as the MERFISH analysis results
% tileNumbers   which tiles we want to put together

% optional input: overlayColours, needs to be a matrix in which the columns define RGB values
if nargin<5
    % set some default colours for the different genes
    overlayColours(1,:) = [0 1 0];      % RGB code for green
    overlayColours(2,:) = [1 0 0];      % RGB code for red
    overlayColours(3,:) = [0 0 1];      % RGB code for blue
    overlayColours(4,:) = [1 1 0];      % RGB code for yellow
    overlayColours(5,:) = [1 0 1];      % RGB code for magenta
    overlayColours(6,:) = [0 1 1];      % RGB code for cyan
    overlayColours(7,:) = [0.9 0.7 0.1];% RGB code for orange
    overlayColours(8,:) = [0.7 1 0.3];  % RGB code for lime green
    overlayColours(9,:) = [0.5 0.5 1];  % RGB code for lavender blue
end

geneLocsPath    = fullfile('D:\MERFISH results',expName,resRound,'geneLocs');          % path for the folder where the XY coordinates for the genes live

%%
% load tiled dapi image
dapiFile = fullfile('D:\MERFISH results',expName,resRound,...
    strcat('tiledDapiImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.tif'));
dapiIm = im2double(imread(dapiFile));

paddedDapi = zeros(size(dapiIm)+2)*NaN;
paddedDapi(2:end-1,2:end-1) = dapiIm;

figure;
imshow(paddedDapi)      % show the dapi in black and white in the background
caxis([0 0.4]);         % change these values to change how bright dapi looks

for g = 1:length(whichGenes)
    hold on                 % this keeps the handle of the figure available for adding to it
    thisColour = cat(3,ones(size(paddedDapi))*overlayColours(g,1),...
        ones(size(paddedDapi))*overlayColours(g,2),...
        ones(size(paddedDapi))*overlayColours(g,3));
    h = imshow(thisColour);
    hold off
    
    thisGene    = whichGenes{g};
    thisGeneFile = fullfile(geneLocsPath,thisGene,...
        strcat(thisGene,'_tiledImage_stacks',num2str(stackNumbers(1)),'to',num2str(stackNumbers(end)),'.tif'));
    thisGeneIm = im2double(imread(thisGeneFile));
    set(h,'AlphaData',thisGeneIm);
    
end

end