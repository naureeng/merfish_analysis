% script that takes merfish TIFF files
% - splits channels
% - creates maximum intensity projections

% based on separateChannels by Yoh Isogai
% circumvents use of ImageJ

% EJ 2019-01-10
% updated 2019-02-14 to rename files

%% initialisations

% startSlice from a blank slate
clear all; close all; clc

% uncomment tic and toc at the end to see how long the script takes to run
% tic

experimentName  = '180712 Merfish';

% set directories
% dataDir     = 'D:\rawData';        % set the path to the folder where the data we want to analyse lives
dataDir     = 'Z:\smfish_data\disk2\Mat';
addpath(dataDir);
thisExpDir  = fullfile(dataDir,experimentName);
ResultsDir  = 'D:\processed';   % path to the folder where the maximum intensity projections in split channels go
if ~exist(fullfile(ResultsDir,experimentName),'dir');
    mkdir(fullfile(ResultsDir,experimentName));
end
thisResultDir = fullfile(ResultsDir,experimentName);

% set parameters
startSlice      = 12;
finishSlice     = 146;
numberSlices    = (finishSlice-startSlice+1)/5;

%% loop through experiment
% get folders of rounds within the experiment
roundFolders    = dir(thisExpDir);                     % this gets a list of all files within the directory
isFolder_rounds = [roundFolders(3:end).isdir];      % in case there are other files in directory, and 1&2 are hidden folders that came up here that we don't care about

% loop to go through round folders
for rr = 1:length(isFolder_rounds)
    
    if isFolder_rounds(rr)
        thisRoundPath   = fullfile(thisExpDir,roundFolders(2+rr).name);
        
        % define bit names here
        switch roundFolders(2+rr).name
            case 'Hyb1' %'Round 1'
                dapiName    = 'hyb1_dapi';
                fidName1    = 'fid_1';
                cy3Name     = 'bit_1';
                cy5Name     = 'bit_5';
                fidName2    = 'fid_5';
                
            case 'Hyb2' %'Round 2'
                dapiName    = 'hyb2_dapi';
                fidName1    = 'fid_2';
                cy3Name     = 'bit_2';
                cy5Name     = 'bit_6';
                fidName2    = 'fid_6';
                
            case 'Hyb3' %'Round 3'
                dapiName    = 'hyb3_dapi';
                fidName1    = 'fid_3';
                cy3Name     = 'bit_3';
                cy5Name     = 'bit_7';
                fidName2    = 'fid_7';
                
            case 'Hyb4' %'Round 4'
                dapiName    = 'hyb4_dapi';
                fidName1    = 'fid_4';
                cy3Name     = 'bit_4';
                cy5Name     = 'bit_8';
                fidName2    = 'fid_8';
                
            case 'Hyb5' %'Round 5'
                dapiName    = 'hyb5_dapi';
                fidName1    = 'fid_9';
                cy3Name     = 'bit_9';
                cy5Name     = 'bit_13';
                fidName2    = 'fid_13';
                
            case 'Hyb6' %'Round 6'
                dapiName    = 'hyb6_dapi';
                fidName1    = 'fid_10';
                cy3Name     = 'bit_10';
                cy5Name     = 'bit_14';
                fidName2    = 'fid_14';
                
            case 'Hyb7' %'Round 7'
                dapiName    = 'hyb7_dapi';
                fidName1    = 'fid_11';
                cy3Name     = 'bit_11';
                cy5Name     = 'bit_15';
                fidName2    = 'fid_15';
                
            case 'Hyb8' %'Round 8'
                dapiName    = 'hyb8_dapi';
                fidName1    = 'fid_12';
                cy3Name     = 'bit_12';
                cy5Name     = 'bit_16';
                fidName2    = 'fid_16';
                
            otherwise   %  then this folder contains something else and we want to skip going through it
                continue  
        end
        
        disp(['loading stacks from ' roundFolders(2+rr).name]);
        
        % get folders of tiles within the round
        tileFolders     = dir(thisRoundPath);
        isFolder_tiles  = [tileFolders(3:end).isdir];
        
        if isempty(find(isFolder_tiles==1))     % then there are no separate folders for the different tiles
            % go directly through the stacks
            dataFiles       = dir(thisRoundPath);
            dataFileNames   = {dataFiles(3:end).name}';
            
            % loop to go through stacks
            for zz = 1:length(dataFileNames)
                %%
                % load files, run maximum intensity projections & save resulting files
                
                % not all files are .tif files, only continue with the ones that are
                if ~isempty(strfind(dataFileNames{zz},'.tif'))
                    
                    thisFilePath = fullfile(thisRoundPath,dataFileNames{zz});
                    info         = imfinfo(thisFilePath);
                    [~,stackName,~] = fileparts(info(1).Filename);
                    if strfind(stackName,'hyb')   % then the stack name also contains the hybridisation round in it, which we want to omit for later file naming purposes
                        stackName = stackName(5:end);
                        if stackName(1) == '_'
                            stackName = stackName(2:end);
                        end
                    end
                    if ~exist(fullfile(ResultsDir,experimentName,stackName,strcat(cy3Name,'.tif')))
                        % in case 
                        disp(['splitting channels & applying maximum intensity projection to ' stackName]);
                        num_images  = numel(info);
                        
                        imageData   = [];
                        images      = {};
                        
                        % maximum intensity projections per channel
                        for j = startSlice:finishSlice
                            images{j} = imread(thisFilePath,j);
                        end
                        
                        % dapi channel
                        imageData.dapi = cat(3,images{startSlice:5:finishSlice});
                        imageData.MAXdapi = max(imageData.dapi,[],3);
                        % fiducial channel - 488
                        imageData.fid = cat(3,images{startSlice+1:5:finishSlice});
                        imageData.MAXfid =  max(imageData.fid,[],3);
                        % cy3 channel - 561
                        imageData.cy3 = cat(3,images{startSlice+2:5:finishSlice});
                        imageData.MAXcy3 =  max(imageData.cy3,[],3);
                        % cy5 channel - 647
                        imageData.cy5 = cat(3,images{startSlice+3:5:finishSlice});
                        imageData.MAXcy5 =  max(imageData.cy5,[],3);
                        
                        % save maximum intensity projections
                        if ~exist(fullfile(thisResultDir,stackName),'dir')
                            mkdir(fullfile(thisResultDir,stackName));
                        end
                        imwrite(imageData.MAXdapi,fullfile(thisResultDir,stackName,strcat(dapiName,'.tif')));
                        imwrite(imageData.MAXfid,fullfile(thisResultDir,stackName,strcat(fidName1,'.tif')));
                        imwrite(imageData.MAXfid,fullfile(thisResultDir,stackName,strcat(fidName2,'.tif')));
                        imwrite(imageData.MAXcy3,fullfile(thisResultDir,stackName,strcat(cy3Name,'.tif')));
                        imwrite(imageData.MAXcy5,fullfile(thisResultDir,stackName,strcat(cy5Name,'.tif')));
                        
                        
                        % to avoid problems, clear the data generated for the next iteration
                        clear thisFilePath info num_images imageData images
                    else
                        disp([stackName ' has already been processed, skipping to next stack']);
                    end
                end
                
            end % end of stack loop
            
            
        else
            % loop to go through tile folders
            warning('Looping through tile folders - this part has not been needed so far, and may be buggy!');
            pause;
            for tt = 1:length(isFolder_tiles)
                thisTilePath    = fullfile(thisRoundPath,tileFolders(2+tt).name);
                
                % get files within the data directory
                dataFiles       = dir(thisTilePath);            % this gets a list of all files within the directory
                dataFileNames   = {dataFiles(3:end).name};      % this extracts their names only into a cell
                
                % loop to go through stacks
                for zz = 1:length(dataFileNames)
                    %%
                    % load files, run maximum intensity projections & save resulting files
                    
                    thisFilePath = fullfile(thisTilePath,dataFileNames{zz});
                    info        = imfinfo(thisFilePath);
                    num_images  = numel(info);
                    
                    imageData   = [];
                    images      = {};
                    
                    % maximum intensity projections per channel
                    for j = startSlice:finishSlice
                        images{j} = imread(thisFilePath,j);
                    end
                    
                    % dapi channel
                    imageData.dapi = cat(3,images{startSlice:5:finishSlice});
                    imageData.MAXdapi = max(imageData.dapi,[],3);
                    % fiducial channel - 488
                    imageData.fid = cat(3,images{startSlice+1:5:finishSlice});
                    imageData.MAXfid =  max(imageData.fid,[],3);
                    % cy3 channel - 561
                    imageData.cy3 = cat(3,images{startSlice+2:5:finishSlice});
                    imageData.MAXcy3 =  max(imageData.cy3,[],3);
                    % cy5 channel - 647
                    imageData.cy5 = cat(3,images{startSlice+3:5:finishSlice});
                    imageData.MAXcy5 =  max(imageData.cy5,[],3);
                    
                    % save maximum intensity projections
                    if ~exist(fullfile(thisResultDir,strcat('Tile_',num2str(tt)),strcat('Stack_',num2str(zz))),'dir')
                        mkdir(fullfile(thisResultDir,strcat('Tile_',num2str(tt)),strcat('Stack_',num2str(zz))));
                    end
                    imwrite(imageData.MAXdapi,fullfile(thisResultDir,strcat('Tile_',num2str(tt)),strcat('Stack_',num2str(zz)),strcat(dapiName,'.tif')));
                    imwrite(imageData.MAXfid,fullfile(thisResultDir,strcat('Tile_',num2str(tt)),strcat('Stack_',num2str(zz)),strcat(fidName1,'.tif')));
                    imwrite(imageData.MAXcy3,fullfile(thisResultDir,strcat('Tile_',num2str(tt)),strcat('Stack_',num2str(zz)),strcat(cy3Name,'.tif')));
                    imwrite(imageData.MAXcy5,fullfile(thisResultDir,strcat('Tile_',num2str(tt)),strcat('Stack_',num2str(zz)),strcat(cy5Name,'.tif')));
                    
                    
                    % to avoid problems, clear the data generated for the next iteration
                    clear thisFilePath info num_images imageData images
                    
                end % end of stack loop
                
            end % end of tile loop
        end
    end % end of "if folder" if loop
    clear dapiName fidName cy3Name cy5Name
end % end of round loop

disp('done.');
% toc
