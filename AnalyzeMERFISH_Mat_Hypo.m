%AnalyzeMERFISH_yohmer.m
%Collection of scripts taken from Zhuang lab matlab codes

%First load imageData and fiducialData
%Then try to come up with fiducial alignments
%Need another script to automate the generation of codebook.

%%
%This section defines parameters.
%

exampleDataPath = '/Users/mathewedwards/Dropbox (UCL - SWC)/merfish_analysis/data/'
dataFolderName = 'merfish_run4_image09'
analysisBasePath = [exampleDataPath,'MERFISH_Demo_Output_run3/']; % This path can be changed to change the location where analysis will be saved

parameters.imageTag = 'STORM'          % Initial name--typically describes the scope
parameters.imageMListType = 'alist'   % Tag for molecule list file for found RNA spots
parameters.fiducialMListType = 'list'  % Tag for molecule list file for found fiducial spots

% Setup parameters of the run
parameters.numHybs = 16                % The number of hybridization/imaging rounds
%parameters.bitOrder = fliplr(1:16);   % The order in which bits are imaged. The example was run in reverse
parameters.bitOrder = 1:16 %canonical order

% Setup parameters for constructing words
parameters.wordConstMethod = 'perLocalization' % A flag to specify the word construction method. This is the preferred method.

% Setup parameters for decoding data
parameters.codebookPath = [exampleDataPath dataFolderName '/merfish_codebook.fasta']           % Insert path to the codebook
parameters.exactMap = CodebookToMap(parameters.codebookPath, ...
        'keyType', 'binStr')
parameters.errCorrFunc = @SECDEDCorrectableWords   % The function used to generate a map for correctable words
%parameters.FPKMData = LoadByteStream(...
%    [exampleDataPath '/FPKM_data/FPKMData.matb'])  % Insert path to the FPKMdata file

% Setup FOV/cells to analyze
parameters.cellsToAnalyze = []       % Analyze all cells if empty

% Setup parameters for saving results
parameters.savePath = SetFigureSavePath(analysisBasePath, 'makeDir', true)

% Configure parameters for generating different reports
parameters.reportsToGenerate = cell(0,2);
parameters.reportsToGenerate(end+1,:) = {'fiducialReport2', 'off'}; % {report name, 'off'/'on' do not/do display figure}
parameters.reportsToGenerate(end+1,:) = {'numOnBitsHistByCell', 'off'}; 
parameters.reportsToGenerate(end+1,:) = {'focusLockReportByCell', 'off'};
parameters.reportsToGenerate(end+1,:) = {'totalFPKMReport', 'off'};
parameters.reportsToGenerate(end+1,:) = {'cellByCellFPKMReport', 'off'};
parameters.reportsToGenerate(end+1,:) = {'cellWithWordsImage', 'off'};
parameters.reportsToGenerate(end+1,:) = {'molStats', 'off'}; 
parameters.reportsToGenerate(end+1,:) = {'molDistStats', 'off'};
parameters.reportsToGenerate(end+1,:) = {'compositeHybImage', 'off'};
parameters.reportsToGenerate(end+1,:) = {'hamming1DReportAllGenes', 'off'};
parameters.reportsToGenerate(end+1,:) = {'bitFlipProbabilitiesAverage', 'off'};
parameters.reportsToGenerate(end+1,:) = {'bitFlipProbabilitiesAllGenes', 'off'};
parameters.reportsToGenerate(end+1,:) = {'confidenceRatioReport', 'off'};

parameters.overwrite = true;                % Overwrite existing files
parameters.figFormats = {'fig', 'png'};     % Output formats
parameters.useSubFolderForCellReport = true; 
parameters.saveAndClose = true; % Save figure once created, then close it

%%
%Load files
%

dataPath = [exampleDataPath dataFolderName]
varargin = {'parameters', parameters}

%% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% Parameters for parsing file names
defaults(end+1,:) = {'imageTag', 'string', 'STORM'};        % Base tag for all images
defaults(end+1,:) = {'imageMListType', 'string', 'alist'};  % Flag for image mlist
defaults(end+1,:) = {'fiducialMListType', 'string','list'}; % Flag for fiducial mlists

% Parameters for parsing file names
defaults(end+1,:) = {'fileExt', 'string', 'bin'};           % Delimiters for bin files
defaults(end+1,:) = {'fieldNames', 'cell', ...              % Labels for fields in image name
    {'movieType', 'hybNum', 'cellNum', 'isFiducial', 'binType'}}; 
defaults(end+1,:) = {'fieldConv', 'cell', ...               % Conversion functions for fields in image name
    {@char, @str2num, @str2num, @(x)strcmp(x,'c2'), @char}};
defaults(end+1,:) = {'appendExtraFields', 'bool', true};    % How to handle names that don't match this pattern

% Parameters for the cells to analyze
defaults(end+1,:) = {'cellsToAnalyze', 'array', []};        % List of cell/FOV ids to analyze

% Parameters on the number of hybridizations
defaults(end+1,:) = {'numHybs', 'nonnegative', 16};         % Number of hybridizations
defaults(end+1,:) = {'bitOrder', 'boolean', 1:16};          % Order of bits

% Parameters for fiducial tracking
defaults(end+1,:) = {'maxD', 'nonnegative', 20};             % Maximum distance for fiducial tracking
defaults(end+1,:) = {'fiducialFrame', 'nonnegative', 1};    % Reference frame for fiducial markers
defaults(end+1,:) = {'fiducialWarp2Hyb1', 'boolean', false};

% Parameters for constructing words from spots
defaults(end+1,:) = {'maxDtoCentroid', 'nonnegative', 5};   % Distance between spots in different rounds

% Parameters for decoding words
defaults(end+1,:) = {'codebookPath', 'filePath', ''};       % Path to codebook
defaults(end+1,:) = {'codebook', 'struct', []};             % Codebook structure
defaults(end+1,:) = {'exactMap', 'map', []};                % containers.Map for decoding exact matches
defaults(end+1,:) = {'correctableMap', 'map', []};          % containers.Map for decoding correctable matches
defaults(end+1,:) = {'errCorrFunc', 'function', ...         % Error correction function
    @SECDEDCorrectableWords};
defaults(end+1,:) = {'keyType', {'int', 'binStr'}, ...      % Display type for binary word, e.g. binary or decimal
    'binStr'};

% Parameters for progress reports and intermediate figures
defaults(end+1,:) = {'savePath', 'fileDir', ''};            % Path to save incidental figures
defaults(end+1,:) = {'reportsToGenerate', 'cell', ...       % List of flags for generating different reports
    cell(0,2)}; 
defaults(end+1,:) = {'verbose', 'boolean', false};          % Display progress?
defaults(end+1,:) = {'printedUpdates', 'boolean', true};    % Display additional forms of progress?

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
%if nargin < 1 || ~(exist(dataPath) == 7) % 7=folder
%    error('matlabFunctions:invalidArguments', 'A valid data path is required.');
%end

%% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename)


%% -------------------------------------------------------------------------
% Import Java utility for generating random IDs
% -------------------------------------------------------------------------
import java.util.UUID;

%% -------------------------------------------------------------------------
% Find data for analysis
% -------------------------------------------------------------------------
if parameters.printedUpdates
    display('--------------------------------------------------------------');
    display(['Finding cells in ' dataPath]);
end

foundFiles = BuildImageDataStructures(dataPath, 'parameters', parameters)
numCells = max([foundFiles.cellNum])

if parameters.printedUpdates
    display(['    Found ' num2str(numCells) ' cells'])
end

if numCells == 0
    error('No valid cells found')
end

% -------------------------------------------------------------------------
% Prepare loop variables
% -------------------------------------------------------------------------
words = [];
totalImageData = [];
totalFiducialData = [];

%% -------------------------------------------------------------------------
% Determine cells to analyze
% -------------------------------------------------------------------------
if isempty(parameters.cellsToAnalyze)
    cellIDs = 1:numCells;
else
    cellIDs = parameters.cellsToAnalyze(parameters.cellsToAnalyze >= 1 & ...
        parameters.cellsToAnalyze <= numCells);
end


% -------------------------------------------------------------------------
% Load codebook and generate maps
% -------------------------------------------------------------------------
if isempty(parameters.codebook) && ~isempty(parameters.codebookPath)
    parameters.codebook = fastaread(parameters.codebookPath)
end

if ~isempty(parameters.codebook)
    parameters.exactMap = CodebookToMap(parameters.codebook, ...
        'keyType', parameters.keyType);
    if ~isempty(parameters.errCorrFunc)
        parameters.correctableMap = CodebookToMap(parameters.codebook, ...
            'keyType', parameters.keyType, ...
            'errCorrFunc', parameters.errCorrFunc, ...
            'mapContents', 'correctable');
    end
end

if parameters.printedUpdates
    display('--------------------------------------------------------------');
    if isempty(parameters.exactMap)
        display('No codebook provided. Found words will not be decoded.');
    elseif isempty(parameters.correctableMap)
        display('No error correction will be applied.');
    end
end

%% -------------------------------------------------------------------------
% Loop over all cells
% -------------------------------------------------------------------------
for i=cellIDs
   % i = 1
    % ---------------------------------------------------------------------
    % Display cell number
    % ---------------------------------------------------------------------
    if parameters.printedUpdates
        display('--------------------------------------------------------------');
        display(['Analyzing data for cell ' num2str(i) ' of ' num2str(numCells)]);
    end
    
    % ---------------------------------------------------------------------
    % Identify all files for this cell
    % ---------------------------------------------------------------------
    imageData = foundFiles(strcmp({foundFiles.movieType}, parameters.imageTag) & ...
        strcmp({foundFiles.binType}, parameters.imageMListType) & ...
        [foundFiles.cellNum] == cellIDs(i) & ...
        ~[foundFiles.isFiducial]);
    
    fiducialData = foundFiles([strcmp({foundFiles.movieType}, parameters.imageTag)] & ...
        strcmp({foundFiles.binType}, parameters.fiducialMListType) & ...
        [foundFiles.cellNum] == cellIDs(i) & ...
        [foundFiles.isFiducial]);

    if parameters.printedUpdates & parameters.verbose
        display(['    Found ' num2str(length(imageData)) ' image files'])
        for j=1:length(imageData)
            display(['       ' imageData(j).filePath])
        end
        display(['    Found ' num2str(length(fiducialData)) ' fiducial files'])
        for j=1:length(imageData)
            display(['       ' fiducialData(j).filePath])
        end
    end
    
    % ---------------------------------------------------------------------
    % Cross checks on found files
    % ---------------------------------------------------------------------
    if length(imageData) ~= length(fiducialData)
        warning('matlabFunctions:AnalyzeMultiFISH', ...
            ['Cell ' num2str(i) ' does not have equal numbers of data and fiducial images']);
        continue;
    end
    if length(imageData) < parameters.numHybs
        warning('matlabFunctions:AnalyzeMultiFISH', ...
            ['Cell ' num2str(i) ' has fewer data images than hybs']);
        continue;
    end
    if length(fiducialData) < parameters.numHybs
        warning('matlabFunctions:AnalyzeMultiFISH', ...
            ['Cell ' num2str(i) ' has fewer fiducial images than hybs']);
        continue;
    end
    
    % ---------------------------------------------------------------------
    % Sort files based on hyb number: Almost certainly not necessary
    % ---------------------------------------------------------------------
    [~, sind] = sort([imageData.hybNum]);
    imageData = imageData(sind);
    [~, sind] = sort([fiducialData.hybNum]);
    fiducialData = fiducialData(sind);
    
    % ---------------------------------------------------------------------
    % Generate and append unique IDs for each image and fiducial data set
    % ---------------------------------------------------------------------
    %for j=1:length(imageData)
    %    imageData(j).uID = char(UUID.randomUUID())
    %end
    %for j=1:length(fiducialData)
    %    fiducialData(j).uID = char(UUID.randomUUID())
    %end
    
    % ---------------------------------------------------------------------
    % Load and transfer information on the corresponding dax: skip for now.
    % ---------------------------------------------------------------------
    %imageData = TransferInfoFileFields(imageData, 'parameters', parameters);
%    fiducialData = TransferInfoFileFields(fiducialData, 'parameters', parameters);
    
    % ---------------------------------------------------------------------
    % Generate a measure of focus lock quality for all images
    % ---------------------------------------------------------------------
    %imageData = GenerateFocusLockQuality(imageData, 'parameters', parameters);

    % ---------------------------------------------------------------------
    % Load Molecule Lists and Fiducial Positions
    % ---------------------------------------------------------------------
    if parameters.printedUpdates & parameters.verbose
        display('--------------------------------------------------------------');
        display('Loading molecules lists');
    end
    for j=1:parameters.numHybs
        imageData(j).mList = ReadMasterMoleculeList(imageData(j).filePath, ...
            'compact', true, 'transpose', true, 'verbose', false);
        fiducialData(j).mList = ReadMasterMoleculeList(fiducialData(j).filePath, ...
            'compact', true, 'transpose', true, 'verbose', false);
        
        if parameters.printedUpdates & parameters.verbose
            display(['    ' imageData(j).name ': ' num2str(length(imageData(j).mList.x)) ' molecules']);
            display(['    ' fiducialData(j).name ': ' num2str(length(fiducialData(j).mList.x)) ' beads']);
        end
    end
     
    % ---------------------------------------------------------------------
    % Create Tiled Image (if desired)
    % ---------------------------------------------------------------------
    %GenerateTiledImage(imageData, 'parameters', parameters);
    
        % ---------------------------------------------------------------------
    % Add Transforms to Fiducial Data
    % ---------------------------------------------------------------------
    [fiducialData, parameters] = AlignFiducials(fiducialData, ...
        'parameters', parameters);
    
    % ---------------------------------------------------------------------
    % Transform Image Data and Transfer Fiducial Data
    % ---------------------------------------------------------------------
    [imageData, parameters] = TransformImageData(imageData,fiducialData,'parameters',parameters);
    
    % -------------------------------------------------------------------------
    % Create Words from Spots
    % -------------------------------------------------------------------------
    [wordsByCell, parameters] = CreateWords(imageData, 'parameters', parameters);

    % -------------------------------------------------------------------------
    % Decode words
    % -------------------------------------------------------------------------
    if ~isempty(parameters.codebook)
       [wordsByCell, parameters] = DecodeWords(wordsByCell, parameters.exactMap, ...
           parameters.correctableMap, 'parameters', parameters);
           if parameters.printedUpdates
                display(['    Found ' num2str(sum([wordsByCell.isExactMatch])) ' exact matches']);
                display(['    Found ' num2str(sum([wordsByCell.isCorrectedMatch])) ' corrected matches']);
           end
    end
    
    % ---------------------------------------------------------------------
    % Create Composite image with words
    % ---------------------------------------------------------------------
    % GenerateCompositeImage(wordsByCell, imageData, 'parameters', parameters);
    
    % ---------------------------------------------------------------------
    % Create Cell By Cell On Bit Histogram
    % ---------------------------------------------------------------------
    % GenerateOnBitHistograms(wordsByCell, 'parameters', parameters);
    
    % -------------------------------------------------------------------------
    % Append Words and imageData
    % -------------------------------------------------------------------------
    words = [words wordsByCell];
    totalImageData = [totalImageData imageData];
    totalFiducialData = [totalFiducialData fiducialData];
end

% change fiducial max search (for images that are significantly shifted)
%parameters.maxD = 100

%% ---------------------------------------------------------------------
    % Add Transforms to Fiducial Data - AlignFiducials worked ok..
    % ---------------------------------------------------------------------
    [fiducialData, parameters] = AlignFiducials(fiducialData, ...
        'parameters', parameters);

%%    ---------------------------------------------------------------------
    % Transform Image Data and Transfer Fiducial Data
    % ---------------------------------------------------------------------
    [imageData, parameters] = TransformImageData(imageData,fiducialData,'parameters',parameters)
    
    % -------------------------------------------------------------------------
    % Create Words from Spots
    % -------------------------------------------------------------------------
    [wordsByCell, parameters] = CreateWords(imageData, 'parameters', parameters)
    
%% % -------------------------------------------------------------------------
    % Decode words
    % -------------------------------------------------------------------------
    if ~isempty(parameters.codebook)
       [wordsByCell, parameters] = DecodeWords(wordsByCell, parameters.exactMap, ...
           parameters.correctableMap, 'parameters', parameters);
           if parameters.zprintedUpdates
                display(['    Found ' num2str(sum([wordsByCell.isExactMatch])) ' exact matches']);
                display(['    Found ' num2str(sum([wordsByCell.isCorrectedMatch])) ' corrected matches']);
           end
    end
    
%% ignore thie part - it won't work unless we have dax file.
    % ---------------------------------------------------------------------
    % Create Composite image with words
    % ---------------------------------------------------------------------
    GenerateCompositeImage(wordsByCell, imageData, 'parameters', parameters);
    
    %% ---------------------------------------------------------------------
    % Create Cell By Cell On Bit Histogram
    % ---------------------------------------------------------------------
    GenerateOnBitHistograms(wordsByCell, 'parameters', parameters);
    
    %% -------------------------------------------------------------------------
    % Append Words and imageData
    % -------------------------------------------------------------------------
    words = [words wordsByCell];
    totalImageData = [totalImageData imageData];
    totalFiducialData = [totalFiducialData fiducialData];
    
    
%%
%Export only the exact matches for this set.
%this set had ~3200 decoded words, ~6400 dots, whereas had 43K correctable
%words, which means that there were something wrong.
    exactMatchWords = [];
    
    for i=1:length(words)
        if words(i).isExactMatch
            exactMatchWords = [exactMatchWords words(i)]
        end
    end
    
  %% This is bad
  %struct2csv(exactMatchWords, '/Users/yohisogai/Desktop/data/merfish180709.csv');
  %exactMatchwordsCell = struct2cell(exactMatchWords);
  %xlswrite('/User/yohisogai/Desktop/data/merfish180711.xls', exactMatchwordsCell);
  %If I want to export the structure array to R, it's best if I save as
  %.mat file and load to R using R.matlab package, readMat(). However, this
  %took a while and it's not recommended. It's best if I focus on what I
  %need to extract the information and load as data.frame. The structure
  %array contains nested lists and it's difficult to handle using R.
  
  %the fiducial images in this set is really screwed up. It's not expected
  %to yield much of anything meaningful. Also, it's possible that I should
  %not set 1 pixel shift from the centre as cutoff. The signal is sparse,
  %and each spot spans 5~6 pixels. Within 1.5 to 2 pixel shift may be more
  %adequate.
  
  
  %Instead I will generate a table by gene, average xc and yc.
  i = 2;
  %Sst
  geneCoordinateList = [];
  wsInds = find(isspace(parameters.codebook(i).Sequence));
  geneName = parameters.codebook(i).Sequence(1:(wsInds(1)-1))

  for j = 1:length(exactMatchWords)
      if strcmp(exactMatchWords(j).geneName, geneName);
       new.name = geneName
       new.xc = nanmean([exactMatchWords(j).xc])
       new.yc = nanmean([exactMatchWords(j).yc])
       geneCoordinateList = [new geneCoordinateList]
      end
  end
  
  figure
  hold on
  scatter([geneCoordinateList.xc], [geneCoordinateList.yc],7,'filled','red')
  xlim([-400 2200])
  ylim([-400 2200])
  pbaspect([1 1 1])
  
  %% second layer
  i = 40;
  %gad1
  colorchoice = 'cyan'
  geneCoordinateList = [];
  wsInds = find(isspace(parameters.codebook(i).Sequence));
  geneName = parameters.codebook(i).Sequence(1:(wsInds(1)-1))

  for j = 1:length(exactMatchWords)
      if strcmp(exactMatchWords(j).geneName, geneName);
       new.name = geneName
       new.xc = nanmean([exactMatchWords(j).xc]);
       new.yc = nanmean([exactMatchWords(j).yc]);
       geneCoordinateList = [new geneCoordinateList];
      end
  end
  scatter([geneCoordinateList.xc], [geneCoordinateList.yc],7, 'filled',colorchoice)
  
  %% third layer
  i = 39;
  %gad2
  colorchoice = 'yellow'
  geneCoordinateList = [];
  wsInds = find(isspace(parameters.codebook(i).Sequence));
  geneName = parameters.codebook(i).Sequence(1:(wsInds(1)-1))

  for j = 1:length(exactMatchWords)
      if strcmp(exactMatchWords(j).geneName, geneName);
       new.name = geneName
       new.xc = nanmean([exactMatchWords(j).xc]);
       new.yc = nanmean([exactMatchWords(j).yc]);
       geneCoordinateList = [new geneCoordinateList];
      end
  end
  scatter([geneCoordinateList.xc], [geneCoordinateList.yc],7, 'filled',colorchoice)
  
  %% fourth layer
  
  i = 7;
  %Oxt
  colorchoice = 'green'
  geneCoordinateList = [];
  wsInds = find(isspace(parameters.codebook(i).Sequence));
  geneName = parameters.codebook(i).Sequence(1:(wsInds(1)-1))

  for j = 1:length(exactMatchWords)
      if strcmp(exactMatchWords(j).geneName, geneName);
       new.name = geneName
       new.xc = nanmean([exactMatchWords(j).xc]);
       new.yc = nanmean([exactMatchWords(j).yc]);
       geneCoordinateList = [new geneCoordinateList];
      end
  end
  scatter([geneCoordinateList.xc], [geneCoordinateList.yc],7, 'filled',colorchoice)
  
    %% fifth layer
  
  i = 33;
  %Bdnf
  colorchoice = 'blue'
  geneCoordinateList = [];
  wsInds = find(isspace(parameters.codebook(i).Sequence));
  geneName = parameters.codebook(i).Sequence(1:(wsInds(1)-1))

  for j = 1:length(exactMatchWords)
      if strcmp(exactMatchWords(j).geneName, geneName);
       new.name = geneName
       new.xc = nanmean([exactMatchWords(j).xc]);
       new.yc = nanmean([exactMatchWords(j).yc]);
       geneCoordinateList = [new geneCoordinateList];
      end
  end
  scatter([geneCoordinateList.xc], [geneCoordinateList.yc],7, 'filled',colorchoice)

  
  %% seventh layer
  
  i = 25;
  %Trh
  colorchoice = 'magenta'
  geneCoordinateList = [];
  wsInds = find(isspace(parameters.codebook(i).Sequence));
  geneName = parameters.codebook(i).Sequence(1:(wsInds(1)-1))

  for j = 1:length(exactMatchWords)
      if strcmp(exactMatchWords(j).geneName, geneName);
       new.name = geneName
       new.xc = nanmean([exactMatchWords(j).xc]);
       new.yc = nanmean([exactMatchWords(j).yc]);
       geneCoordinateList = [new geneCoordinateList];
      end
  end
  scatter([geneCoordinateList.xc], [geneCoordinateList.yc],7, 'filled',colorchoice)
  
  %% export
  set(gca, 'Color', 'none')
  export_fig testhypo2.png -transparent
  