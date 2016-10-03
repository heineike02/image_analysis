%Attempt to use BF to get cell size - not successful
FLUORim = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140524_gluc_drop\Pre\A9-Site_0\img_000000000_RFP_001.tif');
BFoofim = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140524_gluc_drop\Pre\A9-Site_0\img_000000000_BF_000.tif');
BFim = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140524_gluc_drop\Pre\A9-Site_0\img_000000000_BF_001.tif');
BFbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140524_gluc_drop\BG_BF_gluc_drop\img_000000000__000.tif');

BF_BG = BFoofim-BFbg;
BG_BF = BFbg-BFoofim; 

BF_max = max(BF_BG,BG_BF);

[centers,radii,metric] = imfindcircles(BF_max,[5 15], 'ObjectPolarity', 'bright','Sensitivity',.9);...
figure(5); imshow(FLUORim,[]); hold on; viscircles(centers,radii, 'edgecolor', 'b'); title('imfindcircles: ObjPol dark');

figure(1); imagesc(FLUORim); title('FLUORim')
figure(2); imagesc(BFim); title('BFim')
figure(3); imagesc(BFoofim); title('BFoofim')
figure(4); imagesc(BF_max);title('BF_max')



