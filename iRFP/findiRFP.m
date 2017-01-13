function [] = findiRFP()

% load relevant images
imgRFP = imread('t08xy12c1.tif');
imgiRFP = imread('t08xy12c2.tif');
imgComb = imread('t08xy12combined.tif');
imgRFPt0= imread('t01xy12c1.tif');
imgiRFPt0 = imread('t01xy12c2.tif');

% make a good mask
mediRFP=median(double(imgiRFP(:)));
stdiRFP=std(double(imgiRFP(:)));
% binary mask
maskiRFP = imgiRFP>mediRFP+1.5*stdiRFP;
figure(3); imagesc(maskiRFP)
% remove small objects
cc = bwconncomp(maskiRFP);
stats = regionprops(cc, 'Area');
idx=find([stats.Area]>10);
maskiRFP2 = ismember(labelmatrix(cc),idx);
figure(4); imagesc(maskiRFP2);
% fill out holes in objects
maskiRFP3 = imopen(maskiRFP2, strel('disk',1));
figure(5); imagesc(maskiRFP3);
maskiRFP4 = imclose(maskiRFP3, strel('disk',2));
figure(6); imagesc(maskiRFP4);
% remove edge objects
maskiRFP5=imclearborder(maskiRFP4); 
figure(7); imagesc(maskiRFP5);
% find area and centroids and remove small objects
dd = bwconncomp(maskiRFP5);
stats1= regionprops(dd,'Area','Centroid');
idx1=find([stats1.Area]>30 & [stats1.Area]<300);
maskiRFP6 = ismember(labelmatrix(dd),idx1);
figure(8); imagesc(maskiRFP6);
% shift the "nucleus" so it overlaps better with the mCherry localization
ee = bwconncomp(maskiRFP6);
for i = 1:length(ee.PixelIdxList); ee.PixelIdxList{i} = ee.PixelIdxList{i} + 3; end
ee1=labelmatrix(ee); 
figure(9); imagesc(ee1);
maskiRFP7 = ee1>0;
figure(10); imagesc(maskiRFP7);

% apply mask to grayscale image
imgRFP1 = maskiRFP7.*double(imgRFP);
figure(11); imagesc(imgRFP1);
% obtain intensity values for image
ff = bwconncomp(imgRFP1);
stats2= regionprops(ff,imgRFP1, 'Area','Centroid','PixelList', 'PixelValues');

figure(8); imagesc(maskiRFP6);