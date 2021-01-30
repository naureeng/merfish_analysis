% test_cellSegmentationTools

dataFolder = 'C:\Users\Elina\Dropbox (UCL - SWC)\For Elina\oligodT test images';
imName = 'oligodT_p0_m030.tif';                             % randomly load a stack for testing purposes
testIm = imread(fullfile(dataFolder,imName));
dapi = testIm(:,:,1);
oldt = testIm(:,:,2);

%% just see what happens when doing a watershed directly on the raw images
% test1 = watershed(dapi);
% test2 = watershed(oldt);
% 
% % plot results
% figure;
% subplot(2,2,1)
% imagesc(dapi); axis square; axis off; colormap(parula)
% title('dapi')
% subplot(2,2,2)
% imagesc(test1); axis square; axis off; colormap(gray)
% title('watershed dapi')
% 
% subplot(2,2,3)
% imagesc(oldt); axis square; axis off; colormap(parula)
% title('oligo dt')
% subplot(2,2,4)
% imagesc(test2); axis square; axis off; colormap(gray)
% title('watershed oligodt')

% -> this is not very good

%% try converting the dapi into a binary black and white image

h = imhist(dapi);
figure; bar(h)
T = otsuthresh(h);      % uses Otsu thresholding to find the threshold to set for conversion to bw
dapiBW = im2bw(dapi,T);

figure; 
subplot(1,2,1)
imagesc(dapi); axis square; axis off; colormap(parula)
title('dapi');
subplot(1,2,2)
imagesc(dapiBW); axis square; axis off;
title('thresholded BW dapi');

% -> not bad, however cells that are close together are not segmented properly

% try using the edge  detection method
sx = fspecial('sobel');
sy = sx';
gx = imfilter(double(dapi),sx,'replicate');
gy = imfilter(double(dapi),sy,'replicate');
grad = sqrt(gx.*gx + gy.*gy);
grad = grad/max(grad(:));

figure; 
subplot(1,2,1)
imagesc(dapi); axis square; axis off; colormap(parula)
title('dapi');
subplot(1,2,2)
imagesc(grad); axis square; axis off;
title('dapi edges');

hgrad = imhist(grad);
Q = percentile2i(hgrad,0.94);

markerImage = grad > Q;
figure; imshow(markerImage)

dapip = double(dapi).*double(markerImage);
figure; imagesc(dapip); axis square; axis off; colormap(parula)

hp = imhist(uint8(dapip));
hp(1) = 0;
figure; bar(hp);
Tp = otsuthresh(hp);
g = im2bw(dapi, Tp);
figure; 
subplot(1,3,1); imagesc(dapi); axis square; axis off;
subplot(1,3,2); imagesc(dapiBW); axis square; axis off;
subplot(1,3,3); imagesc(g); axis square; axis off;

% -> better

%% watershed

dbw = im2bw(dapi, graythresh(dapi));
dbw_c = ~dbw;
D = bwdist(dbw_c);
L = watershed(-D);
w = L == 0;

% to compare the different thresholded images
% figure;
% subplot(2,3,1); imagesc(dapi); axis square; axis off; colormap(gray)
% subplot(2,3,4); imagesc(dapiBW); axis square; axis off;
% subplot(2,3,5); imagesc(g); axis square; axis off;
% subplot(2,3,6); imagesc(dbw); axis square; axis off;

figure; 
subplot(2,2,1); imagesc(dapi); axis off, axis square
subplot(2,2,2); imagesc(dbw); axis off, axis square
subplot(2,2,3); imagesc(D); axis off, axis square
subplot(2,2,4); imagesc(w); axis off, axis square

segmDapi = dbw & ~w;
figure; 
subplot(1,2,1); imagesc(dapi), axis off, axis square, title('dapi')
subplot(1,2,2); imagesc(segmDapi), axis off, axis square, title('segmented dapi')


% try all the different threshold images in the watershed
figure;
subplot(2,3,1); imagesc(dapi); axis off, axis square, colormap(gray)
title('raw dapi');
for tt = 1:3
    switch tt
        case 1
            thisTdapi = dbw;
            ttl = 'graythreshold';
        case 2
            thisTdapi = dapiBW;
            ttl = 'simple otsuthreshold';
        case 3
            thisTdapi = g;
            ttl = 'otsu threshold on borders';
    end
    Tdapi_c = ~thisTdapi;
    D = bwdist(Tdapi_c);
    L = watershed(-D);
    w = L == 0;
    segmDapi = thisTdapi & ~w;
    subplot(2,3,3+tt); imagesc(segmDapi); axis off, axis square, title(ttl)
end

% -> graythresh and simple otsu equally as good

% get boundaries

[boundaries, labels, n] = bwboundaries(segmDapi,'noholes');
nPixPerROI = cellfun('length',boundaries);                  % this gives how many pixels are within boundarie of each area
errors = find(nPixPerROI<50 | nPixPerROI > 150);
test = labels;
for e = 1:length(errors)
    thisVal = errors(e);
    test(test==thisVal) = 0;
end
% -> not optimal, what we really want is the number of pixels inside boundaries

errors = [];
cE = 0;
for ii = 1:n
    thisPutCell = boundaries{ii};
    area = polyarea(thisPutCell(:,2),thisPutCell(:,1));
    if area < 200 | area > 1500
        cE = cE+1;
        errors(cE) = ii;
    end
end
test = labels;
for e = 1:length(errors)
    thisVal = errors(e);
    test(test==thisVal) = 0;
end

% % add internal & external markers
% dapiO = imcomplement(dapi);
% mins = imregionalmin(dapiO);


%%
% nbS = 5;
% dapiSTDfilt = stdfilt(dapi, ones(nbS));
% % dapiSTDgaussFilt = imgaussfilt(dapiSTDfilt,1);
% figure; 
% subplot(1,2,1)
% imagesc(dapi); axis square; axis off; colormap(parula)
% title('dapi');
% subplot(1,2,2)
% imagesc(dapiSTDfilt); axis square; axis off;
% title(['dapi STD (neighbourhood size ' num2str(nbS) ')']);
% % subplot(1,3,3)
% % imagesc(dapiSTDgaussFilt); axis square; axis off;
% % title(['dapi STD with gauss filt']);
% % -> adding the gaussian filter doesn't really help
% % h2 = imhist(dapiSTDfilt);
% % 
% % g = localthresh(oldt, ones(nbS), 10, 1, 'global'); figure; imagesc(g)
% % -> doesn't work either
