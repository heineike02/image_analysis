function [STATS] = bfanalyze(BFim, FLUORim, shift);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION - Does cell segmentation for a single image                   %
%                                                                       %
%INPUTS:    BF - brightfield image name                                 %
%           FLUOR - fluorescence image name                             %
%                                                                       %
%OUTPUTS:   STATS - metrics for each object identified in image         %
%                   'Area' - area of each cell found                    %
%                   'Centroid' - XY coordinates of cell                 %
%                   'PixelIdxList' - linear indices of pixels           %
%                   'PixelList' - locations of pixels in each cell      %
%                   'MeanIntensity' - mean intensity of each cell       %
%                   'PixelValues' - intensities of each pixel in cell   %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find circles strategy
%im = imread('glucgal2T0XY01C2.tif');
[centers,radii,metric] = imfindcircles(BFim,[3 13], 'ObjectPolarity', 'dark');...
figure; imshow(BFim,[]); hold on; viscircles(centers,radii, 'edgecolor', 'b'); title('imfindcircles: ObjPol dark');

%align the fluorescence to brightfield images
%imf = imread('glucgal2T02XY01C2.tif');
if shift == 1;
    centersx = centers(1:end,1)-6; 
    centersy = centers(1:end,2)+10;
    centersnew= [centersx,centersy];
else
    centersnew = centers;
end

%figure; imshow(FLUORim,[]); hold on; viscircles(centersnew,radii, 'edgecolor', 'b'); title('imfindcircles: ObjPol dark');


%store circles to compress into a binary mask
masktemp=zeros(1952,1952,length(centersnew(1:end,1)));
for i = 1:length(centersnew(1:end,1))
    BWbg = zeros(1952, 1952);
    %generate a circle mask
    radius =  radii(i);
    %circ_mask = double(getnhood(strel('ball',radius,radius,0)));
    %figure; imshow(circ_mask, 'InitialMagnification', 500)
    BWbg(round(centersnew(i,2)), round(centersnew(i,1))) = 1;
    masktemp(:,:,i)= conv2(BWbg, fspecial('disk',radii(i)), 'same');
    %masktemp= conv2(BWbg, fspecial('disk',radii(i)), 'same'); 
end

%make a binary mask
mask_bw = sum(masktemp,3) >0; %figure; imshow(mask_bw,[]);

%multiply binary mask with fluorescence image
imf_masked = double(FLUORim).*mask_bw; %imshow(imf_masked,[]);

%get region properties
L = bwlabel(imf_masked,8); figure; imshow(L,[]);
STATS = regionprops(L, imf_masked, 'Area', 'Centroid', 'PixelIdxList','PixelList', 'MeanIntensity', 'PixelValues');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STUFF DIDN'T USE
% %imread 
% %im = imread();
% figure; imshow(im,[]);
% 
% %filter image
% H= fspecial('disk',2); IMFIL = imfilter(im,H);
% figure; imshow(IMFIL,[]);
% 
% %find edges
% BWed= edge(IMFIL,'sobel'); figure; imshow(BWed,[]);
% 
% %dilate image
% SE=strel('disk',5); BWdil= imdilate(BWed,SE); figure; imshow(im.*BWdil,[]);
% 
% %fill in the small black areas
% SE=strel('disk',5); BWclo= imclose(BWdil,SE); figure; imshow(im.*BWclo,[]);
% 
% %erode the image
% SE=strel('disk',5); BWero= imerode(BWclo,SE); figure; imshow(im.*BWero,[]);