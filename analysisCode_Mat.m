% create variables

expName     = '180712 Merfish';
resRound    = '001';

% to get a list of all genes
geneFolders = dir('D:\MERFISH results\180712 Merfish\geneLocs');
geneList = {geneFolders(3:end).name};


stackNumbers = 2:26;

tileSize = [5 5];
tilingOption = 'snake';

[tileOrder] = createTileOrder(tileSize,tilingOption);

createTiledImage(expName,resRound,geneList,stackNumbers,tileOrder);

%%

whichGenes = {};
whichGenes{1} = 'OXT';
whichGenes{2} = 'SST';

multiGeneOverlays(expName,whichGenes,stackNumbers);