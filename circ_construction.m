%circ - old circ.mat from Jacob
%just saving this as the 1x SC Circ

load 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\BMH_Image_Analysis\circ.mat'
fname_out = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\BMH_Image_Analysis\circSC_1x.mat' ;
save(fname_out, 'circ'); 

%To Do: 
%1)Either figure out how Jacob made the old circ (where did the range come
%from) or just test out our circles with the range we have from these images. 
%2)Make a new 1x Circ for S.C. to replace circ
%3)The border around these are wierd - maybe use a 

%circKL - taken from an RFP image of ?what strain? ?what exposure? on
%16JUL14. 1x magnification.
%fname_in = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\WetLab\20130716\Gluc_pre_ind\Pos5\img_000000000_RFP_000.tif' ;
%fname_out = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\BMH_Image_Analysis\circKL_1x.mat' ;
%y_coord = 198:212;
%x_coord = 80:94;


%circKL_1x - taken from an RFP image of ?what strain? ?what exposure? on
%16JUL14. 1x magnification.
%fname_in = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\WetLab\20130716\Gluc_pre_ind\Pos5\img_000000000_RFP_000.tif' ;
%fname_out = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\BMH_Image_Analysis\circKL_1x.mat' ;
%y_coord = 198:212;
%x_coord = 80:94;

%circSC_15x - taken from an RFP image of SC MSN2(dzf)-mCherry (HES 11-38)
%with 100ms exposure on intensity 8.  
%fname_in = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140204\Gluc_Rep_Post_2\RFP_p3_t11.tiff' ;
%fname_out = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\BMH_Image_Analysis\circSC_15x.mat' ;
%y_coord = 46:71;
%x_coord = 266:291;

%circKL_15x - taken from an YFP image of KL MSN2(dzf)-mCherry (made with plasmid BMH 028)
%with 100ms exposure on intensity 8.  
fname_in = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140204\Gluc_Drop_Pre\YFP_p6_t2.tiff' ;
fname_out = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\BMH_Image_Analysis\circKL_15x.mat' ;
y_coord = 125:143;
x_coord = 265:283;

im = imread(fname_in);

load 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\BMH_Image_Analysis\circ.mat'
circ_out = double(im(y_coord,x_coord));
local_smooth = 3; 
circ_out = conv2(circ_out,ones(local_smooth),'same');
%make circ_out have the same range as circ.  Where did that range come from??

circ_out = circ_out/512^2;
circ_out = circ_out*( max(max(circ)) -  min(min(circ)))/( max(max(circ_out)) -  min(min(circ_out)));
circ_out = circ_out-min(min(circ_out))+min(min(circ));
circ = circ_out
save(fname_out, 'circ'); 