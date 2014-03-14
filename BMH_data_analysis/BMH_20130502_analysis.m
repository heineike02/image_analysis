di = '/Users/ucsf/Google Drive/UCSF/ElSamad_Lab/PKA/WetLab/20130502/Exp2_1_1/Pos0'
storim = 1
mchan = {'RFP', 'YFP'}

bf_mask = 'img_000000000_Default_000.tif'
metadata = 'metadata.txt'


S=load('/Users/ucsf/Google Drive/HES_Lab_Share/Scripts/JSO_Image_Analysis/circ.mat')
circ=S.circ;
circ=circ(1:13,2:14);

im2 = imadjust(im);
level = graythresh(im2)
im3 = im2bw(im2,level)
im4 = medfilt2(im2)
SE = strel('disk',3);
im5 = imopen(im4,SE);
im6 = imdilate(im5,SE); 
im7 = uint16(im6).*im4;
im8 = watershed(1-im7);
im9 = imread([di,'/img_000000000_Default_000.tif']);
im10=deconvlucy((im4),circ,10);

%ImageJ
median filter
threshold (very arbitrary)
Open
Close
Watershed
Analyze particles 
size between 30-300

%what next?  