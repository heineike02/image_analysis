%Note: works with image_analysis v2.0 (previously CovMapCellsBMH20130716)
path(path,'..\..\image_analysis')

% Incorporate Cell Tracking
% Try with K. Lactis
% Visualize individual and mean nuclear localization for glucose shift
% before and after Glucose deprivation in S.C. and K. Lactis
% Repeat for BPAC (in another file)

%load image - first image is S. Cerevisiae pre-induction
%Windows: 
imdir = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20130716\Gluc_Starvation_1\Pos5\'
%Mac
%imdir = '/Users/ucsf/Google Drive/UCSF/ElSamad_Lab/PKA/WetLab/20130716/Gluc_pre_ind/Pos1/'
im = imread([imdir,'img_000000000_RFP_000.tif']);
bf_mask = imread([imdir,'img_000000000_BF_001.tif']);
storeim = 1;
%imbg = double(imread('..\RFP_BG.tif'));
imbg = 1;
load ..\circ.mat
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
%0.16 seems to work best for k.lactis
std_thresh = 0.16;

max_shift_rad = 4;
min_dist_thresh = 3;

for nn = [1]

    pos = num2str(nn)+4;

    %Earlier time: 
    imdir_pre = ['C:\Users\Ben\Box Sync\Data\PKA_Project\20130716\Gluc_Starvation_1\Pos',pos,'\'];
    times_pre = [-120-2*110,-120-110];
    images_pre = {'img_000000000_RFP_000.tif','img_000000001_RFP_000.tif'};
    images_bf_pre = {'img_000000000_BF_001.tif','img_000000001_BF_001.tif'};


    timecoursedata_pre = CovMapCellsBMH(imdir_pre ,images_pre, {'RFP'}, circ, images_bf_pre , imbg , siz , storeim, rad, std_thresh, times_pre);
    timecoursedata_out_pre = MinDist_BMH(imdir_pre, timecoursedata_pre, max_shift_rad, min_dist_thresh, storeim);

     
    x_vec_pre = zeros(length(timecoursedata_out_pre(1).celldata),length(timecoursedata_out_pre));
    y_vec_pre = x_vec_pre;
    nf_vec_pre = x_vec_pre;
    for jj = 1:length(timecoursedata_out_pre);
        x_vec_pre(:,jj) = [timecoursedata_out_pre(jj).celldata.Cxloc];
        y_vec_pre(:,jj) = [timecoursedata_out_pre(jj).celldata.Cyloc];
        nf_vec_pre(:,jj) = [timecoursedata_out_pre(jj).celldata.nf];
    end


    imdir = ['C:\Users\Ben\Box Sync\Data\PKA_Project\20130716\Gluc_Starvation_2\Pos',pos,'\'];
    D = dir(imdir);
    kk = 0;
    kk_bf = 0;
    for jj = 1:length(D)
        if strfind(D(jj).name,'RFP') == 15
            kk = kk+1;
            images{kk} = D(jj).name;
        end
        if strfind(D(jj).name,'BF') == 15
            kk_bf = kk_bf+1;
            images_bf{kk_bf} = D(jj).name;
        end
    end

    %it was about 1:15 between images for this data
    times = [0:75:59*75];

    timecoursedata = CovMapCellsBMH(imdir ,images, {'RFP'}, circ, images_bf , imbg , siz , storeim, rad, std_thresh, times);
    timecoursedata_out = MinDist_BMH(imdir, timecoursedata ,max_shift_rad, min_dist_thresh, storeim);



    x_vec = zeros(length(timecoursedata_out(1).celldata),length(timecoursedata_out));
    y_vec = x_vec;
    nf_vec = x_vec;
    for jj = 1:length(timecoursedata_out);
        x_vec(:,jj) = [timecoursedata_out(jj).celldata.Cxloc];
        y_vec(:,jj) = [timecoursedata_out(jj).celldata.Cyloc];
        nf_vec(:,jj) = [timecoursedata_out(jj).celldata.nf];
    end






    figure(1)
    subplot(1,2,nn)
    plot(x_vec',y_vec')
    title('Position of cells')

    figure(2)
    subplot(1,2,nn)
    clf 
    hold on
    plot(times,nf_vec)
    plot(times_pre,nf_vec_pre)
    xlabel('time (s)')
    ylabel('Nuclear Localization')
    title('Nuclear localization - Unfiltered')

    figure(3)
    subplot(1,2,nn)
    clf
    hold on
    nan_ind = ~isnan(x_vec(:,end));
    nf_vec_filt = nf_vec(nan_ind,:);
    plot(times,nf_vec_filt);

    nan_ind_pre = ~isnan(x_vec_pre(:,end));
    nf_vec_filt_pre = nf_vec_pre(nan_ind_pre,:);
    plot(times_pre,nf_vec_filt_pre)
    title('Nuclear localization - filtered')
    xlabel('time (s)')
    ylabel('Nuclear Localization')

    figure(4)
    subplot(1,2,nn)
    clf
    hold on
    errorbar(times,mean(nf_vec_filt),std(nf_vec_filt))
    errorbar(times_pre,mean(nf_vec_filt_pre),std(nf_vec_filt_pre))
    title('Mean Nuclear Localization')
    xlabel('time (s)')
    ylabel('Nuclear Localization')
    
    nf_vec_mean(:,nn) = mean(nf_vec_filt);
      
end





%bf_mask = 'img_000000000_Default_000.tif';
%metadata = 'metadata.txt'

%S=load('/Users/ucsf/Google Drive/HES_Lab_Share/Scripts/JSO_Image_Analysis/circ.mat')
%circ=S.circ;
%circ=circ(1:13,2:14);
