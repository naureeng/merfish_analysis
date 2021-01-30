function extractGeneCoordinates_perTile(words,geneNameList,exactMatchesOnly)

% function that takes the words structure that is the result of the MERFISH
% analysis & creates a table with gene location coordinates within the
% given frame that the words results comes from for the genes contained
% within geneNameList

% inputs:
% geneNameList      list of the genes we want to know the corrected x y
%                   coorinates for
% words             structure that comes from the MERFISh analysis

%%
% set it as a a default to only use exact matches
if nargin < 3
    exactMatchesOnly = true;
end

if nargin < 2
    disp('no gene names provided');
    addGene = logical(input('Add gene name? 1 for yes '));
    cc = 0;
    while addGene
        cc = cc+1;
        geneNameList{cc} = input('Provide gene name (as it appears in the MERFISH words results structure) ');
        addGene = logical(input('Add another gene? 1 for yes '));
    end
    clear cc
end

nGenes = length(geneNameList);
numHyb = words.numHyb;

pathParts = strsplit(words(1).imagePaths{1},'\');
expName   = char(pathParts{3});
stackName = char(pathParts{4});
resFolName= char(pathParts{5});

pixelDimensions = [2048 2048];  % size of the original tif files in pixels
pixelSize       = 0.105;        % size of a pixel in um

%%
% loop through genes & extract the coordinates, then save
for g = 1:nGenes
   thisGene = char(geneNameList(g));
   msg = ['extracting XY coordinates for ' thisGene ' in ' stackName];
   disp(msg);
   % find indices within words list for which correspond to the gene of interest
   geneIndsC = strfind({words.geneName},thisGene);
   geneInds  = find(not(cellfun('isempty',geneIndsC)))';
   clear geneIndsC
   
   exactMatchInds = find([words.isExactMatch])';
   exactGeneInds = intersect(exactMatchInds,geneInds);
   nSpots        = length(exactGeneInds);
   clear exactMatchInds
   
   roundCoordsX = reshape([words(exactGeneInds).xc],numHyb,length(exactGeneInds))';
   roundCoordsY = reshape([words(exactGeneInds).yc],numHyb,length(exactGeneInds))';
   % sometimes, some spots get a spurious location assigned; negative or
   % outside the image range - exclude those
   for n = 1:nSpots
       if nanmean(roundCoordsX(n,:)) < 1 || nanmean(roundCoordsY(n,:)) < 1
           warning('some locations are negative, setting these to NaN...');
           roundCoordsX(n,:) = NaN;
           roundCoordsY(n,:) = NaN;
       end
       if nanmean(roundCoordsX(n,:)) > pixelDimensions(1) || nanmean(roundCoordsY(n,:)) > pixelDimensions(1)
           warning('some locations are bigger than the pixel Dimensions, setting these to NaN...');
           roundCoordsX(n,:) = NaN;
           roundCoordsY(n,:) = NaN;
       end
   end
   
   Xcoords = nanmean(roundCoordsX,2);
   Xcoords = Xcoords(isfinite(Xcoords));
   Ycoords = nanmean(roundCoordsY,2);
   Ycoords = Ycoords(isfinite(Ycoords));
   
   stackIm = zeros(pixelDimensions);    % initialise the 2D picture of the gene locations we're generating
   
   % to turn the XY coordinates into a 2D matrix, we need integers
   % as a first approximation, round the coordinates to integers - however this will lead to a loss of resolution
   intXc = round(Xcoords);
   intYc = round(Ycoords);
   stackIm(sub2ind(size(stackIm),intXc,intYc)) = 1;         % creates a matrix where at the specified X&Y positions the value is 1
   
   % uncomment if want to figure out whether this shift is more of less than the std of the estimated position
%    % figure out in um how much we shifted the estimated position by rounding, 
%    % and whether this is more or less than the std of the estimated positions between rounds
%    Xshift = [Xcoords - intXc];    % gives shift in um
%    Yshift = [Ycoords - intYc];
%    Xstd = nanstd(reshape([words(exactGeneInds).xc],numHyb,length(exactGeneInds)))';
%    Ystd = nanstd(reshape([words(exactGeneInds).yc],numHyb,length(exactGeneInds)))';
   
   filename = char(strcat(stackName,'_',thisGene,'_XYcoordinates'));
   thisDir = char(fullfile('D:','MERFISH results',expName,resFolName,'geneLocs',thisGene));
   if ~exist(thisDir,'dir')
       mkdir(thisDir);
   end
   save(fullfile(thisDir,filename),'Xcoords','Ycoords','stackIm');

end

end