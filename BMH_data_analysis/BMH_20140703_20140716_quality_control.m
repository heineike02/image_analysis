ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\';
path(ipdir,path)
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20140703_39y_GD_sorb_p5M\';
phase = 'Post\';
site_vec = 1:3;
time_vec = 0:14;
clear celldata_full
celldata_full.RFP = cell(length(site_vec),length(time_vec));
celldata_full.YFP = cell(length(site_vec),length(time_vec));
for ss = 1:length(site_vec)
    for tt = 1:length(time_vec)
        location = ['A11-Site_',int2str(site_vec(ss)),'\'];
        time = sprintf('%02.2d', time_vec(tt));
        channels = {'RFP','YFP'};

        for jj = 1:2
            channel = channels{jj};
            imname = ['img_0000000',time,'_',channel,'_001.tif'];
            im.(channel) = imread([base_dir,phase,location,imname]);
            %load circ
            load([ipdir, 'circKL_15x.mat'])

            bf_mask = []
            imbg_temp = imread([base_dir,'\BG\img_000000000_',channel,'_001.tif']);
            imbg_temp = double(imbg_temp);
            coarse_smooth = 25;
            imbg.(channel) = medfilt2(imbg_temp,[coarse_smooth,coarse_smooth],'symmetric');
            siz = [18,18];
            storeim = 1; 
            rad = 8;
            std_thresh = 0.16; 

            celldata.(channel) = FindCellsBMH(im.(channel), circ, bf_mask , imbg.(channel) , siz, storeim, rad, std_thresh);

            figure(jj)
            clf
            hold on
            imagesc(im.(channel))
            plot([celldata.(channel).Cyloc],[celldata.(channel).Cxloc],'rx','MarkerSize',5,'LineWidth',3')
            title(channel)

        end

    %try with L = rad
    L = rad;
    [celldata_RFP_matched,celldata_YFP_matched] = match_two_channels(celldata.RFP,celldata.YFP,L);
    
    figure(3)
    plot([celldata_RFP_matched.nf],[celldata_YFP_matched.nf],'x')
    
    ss
    tt
    
    
    celldata_full(ss,tt).RFP = celldata_RFP_matched;
    celldata_full(ss,tt).YFP = celldata_YFP_matched;

    end
end

figure(4)
clf 
hold on
for tt = 1:14
    plot([celldata_full(tt).RFP.nf],[celldata_full(tt).YFP.nf],'x')
end

%do this for a time series of images.  Set location to average location for
%each image. track cells and compare nuclear localization in single cells. 


return

%%Need to add the position as a piece of data
fname_save = '20140703_processed_data_KL_RFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec');
tracks_RFP = all_tracks_vec{1};

fname_save = '20140703_processed_data_KL_YFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec');

tracks_YFP = all_tracks_vec{1};

%For a given time point, plot the points for RFP and YFP on the correct
%image

time = 1;
phase = 'Post';


