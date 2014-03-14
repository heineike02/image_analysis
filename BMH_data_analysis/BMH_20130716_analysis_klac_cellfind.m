%Note: works with CovMapCellsBMH20130716

% Incorporate Cell Tracking
% Try with K. Lactis
% Visualize individual and mean nuclear localization for glucose shift
% before and after Glucose deprivation in S.C. and K. Lactis
% Repeat for BPAC (in another file)

%load image - first image is S. Cerevisiae pre-induction
%Windows: 
imdir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\WetLab\20130716\Gluc_pre_ind\Pos5\'
%Mac
%imdir = '/Users/ucsf/Google Drive/UCSF/ElSamad_Lab/PKA/WetLab/20130716/Gluc_pre_ind/Pos1/'
im = imread([imdir,'img_000000000_RFP_000.tif']);
bf_mask = imread([imdir,'img_000000000_BF_001.tif']);
storeim = 1;
%imbg = double(imread('..\RFP_BG.tif'));
imbg = 1;
%load ..\circ.mat
circKL = double(im(198:212,80:94));
local_smooth = 3; 
circKL=conv2(circKL,ones(local_smooth),'same');
%make circKL have the same range as circ:
circKL = circKL/512^2;
circKL = circKL*( max(max(circ)) -  min(min(circ)))/( max(max(circKL)) -  min(min(circKL)));
circKL = circKL-min(min(circKL))+min(min(circ));
save('circKL.mat', 'circKL'); 

%[celldata,N] = FindCellsBMH(im, circ, bf_mask, im_bg , siz, storeim)
%Size of image of individual cell to analyze
siz = [15,15];
%radius of cell for nuclear enrichment calculation
rad = 4;
%std thresh for calling a peak to find a cell
std_thresh = 0.2;

max_shift_rad = 4;
min_dist_thresh = 3;

celldata = FindCellsBMH(im, circKL, bf_mask,imbg, siz, storeim, rad, std_thresh);

nloc = [celldata.nf];
nmed = [celldata.nmi];
xloc = [celldata.Cxloc];
yloc = [celldata.Cyloc];


figure(2)
hist(nloc)

figure(3)
hist(nmed)

celldata1 = FindCellsBMH(im, circ, bf_mask,imbg, siz, storeim, rad, 0.16);

nloc1 = [celldata1.nf];
nmed1 = [celldata1.nmi];
xloc1 = [celldata1.Cxloc];
yloc1 = [celldata1.Cyloc];

figure(4)
hist(nloc1)

figure(5)
hist(nmed1)

figure(1)
imagesc(im)
hold on
pause
plot(yloc,xloc,'xy')
pause
plot(yloc1,xloc1,'ok')

