function [outMask, imgRFP, imgiRFP] = makeMask(maskImg, intImg, stdConst, remS1, remS2, disc, mvDown, showFig)

% load images
imgRFP = intImg;
imgiRFP = maskImg;
% make a good mask
mediRFP=median(double(imgiRFP(:)));
stdiRFP=std(double(imgiRFP(:)));
% binary mask
maskiRFP = imgiRFP>mediRFP+stdConst*stdiRFP;
% remove small objects
cc = bwconncomp(maskiRFP);
stats = regionprops(cc, 'Area');
idx=find([stats.Area]>remS1);
maskiRFP2 = ismember(labelmatrix(cc),idx);
% fill out holes in objects
maskiRFP3 = imopen(maskiRFP2, strel('disk',disc(1)));
maskiRFP4 = imclose(maskiRFP3, strel('disk',disc(2)));
% remove edge objects
maskiRFP5=imclearborder(maskiRFP4); 
% find area and centroids and remove small objects
dd = bwconncomp(maskiRFP5);
stats1= regionprops(dd,'Area','Centroid');
idx1=find([stats1.Area]>remS2(1) & [stats1.Area]<remS2(2));
maskiRFP6 = ismember(labelmatrix(dd),idx1);
% shift the "nucleus" so it overlaps better with the mCherry localization
ee = bwconncomp(maskiRFP6);
for i = 1:length(ee.PixelIdxList); ee.PixelIdxList{i} = ee.PixelIdxList{i} + mvDown; end
ee1=labelmatrix(ee); 
maskiRFP7 = ee1>0;
% figures
if showFig
    figure(1); imagesc(imgRFP);
    figure(2); imagesc(imgiRFP);
    figure(3); imagesc(maskiRFP);
    %figure(4); imagesc(maskiRFP2);
    %figure(5); imagesc(maskiRFP3);
    %figure(6); imagesc(maskiRFP4);
    %figure(7); imagesc(maskiRFP5);
    %figure(8); imagesc(maskiRFP6);
    %figure(9); imagesc(ee1);
    %figure(10); imagesc(maskiRFP7);
end

outMask = maskiRFP7;
end