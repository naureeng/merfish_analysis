% script that loads multilayers tif file into matlab & displays it
% to allow for inspection

geneLocsPath    = fullfile('D:\MERFISH results',expName,'geneLocs');          % path for the folder where the XY coordinates for the genes live

thisGene    = 'AGRP';
stackNums   = 2:26;

imName = fullfile(geneLocsPath,thisGene,...
    strcat(thisGene,'_tiledImage_stacks',num2str(stackNums(1)),'to',num2str(stackNums(end)),'.tif'));

I = imread(imName);

imshow(I);