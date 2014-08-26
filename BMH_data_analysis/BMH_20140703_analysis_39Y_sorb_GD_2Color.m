function all_tracks = BMH_20140703_analysis_39Y_sorb_GD_2Color()

profile on
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

% Filename convention
% 'JSO','Micromanager', or 'HCS_Nikon'
% For example of Micromanager filename convention see BMH_20140127 analysis
% files
fname_conv = 'Micromanager'

%imdirPhase.Pre = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_pre_1\'
%imdirPhase.Post = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_post_1\'

base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20140703_39y_GD_sorb_p5M\'
imdirPhase.Pre = [base_dir,'Pre\']
imdirPhase.Post = [base_dir,'Post\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%input information from experiment here
species = 'KL' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species
species_cell = {'SC','KL'} 

channels = {'BF','RFP','YFP'}
channels_to_image = {'RFP','YFP'}

fname_saveSP.KL = '20140703_processed_data_KL_2color.mat';
%fname_saveSP.SC = '20140703_processed_data_SC.mat'; 

phases =  {'Pre','Post'} %,'Post'} 
shift_timing = [0,10.5]    
%These are the absolute times at which each phase starts.
%timestep method
%time_calc:  Tells the program how to calculate each time value.  If the
%input is a number then that is interpreted as a dt between images
%If the input is a filename, then that is interpreted as metadata that gives
%exact time values.  

%time_calc_phase.Pre = 5
%time_calc_phase.Post = 4
%time_calc_phase.Pre = [imdirPhase.Pre, 'acqdat.txt']
%time_calc_phase.Post = [imdirPhase.Post, 'acqdat.txt']
%extract actual times from metadata in micromanager images
%addpath('C:/Users/Ben/Documents/GitHub/image_analysis/jsonlab');
%time_calc_phase.Pre  = 'metadata.txt'
%time_calc_phase.Post = 'metadata.txt'
generate_metadata_parsed = 0;
metadata_conv_fname = 'C:\Users\Ben\Documents\GitHub\image_analysis\times_from_umanager_metadata.py';
if generate_metadata_parsed ==1;
    for ph = [1:length(phases)]
          phase = phases{ph};
          imdir = imdirPhase.(phase);
          system(['python ', metadata_conv_fname, ' ', imdir])
    end
end
time_calc_phase.Pre =  'metadata_parsed.txt';
time_calc_phase.Post = 'metadata_parsed.txt';


%std thresh for calling a peak to find a cell
%0.16 seemed to work best for k.lactis on 16JUL images
std_threshSP.KL = 0.16;
std_threshSP.SC = 0.2;

%Tracking parameters
%estimated maximum number of tracks
maxdisp_1x = 4;

% optical amplification 
% 1x or 1.5x
op_amp =  '1.5x'
storeim = 1;

bgimg = 0; 
%1 if you have a background image, 0 if you don't. 
%Gets background image depending on channel to image
%Collect background images using micromanager 
if bgimg == 0
    imbg = 1; % default
else
    %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
    if strcmp(channel_to_image,'RFP')
         imbg = imread([base_dir,'\BG\img_000000000_RFP_001.tif']);
    elseif strcmp(channel_to_image,'YFP')
         imbg = imread([base_dir,'\BG\img_000000000_YFP_001.tif']);
         %imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\YFP_BG\img_000000000_Default_000.tiff');
    else
         'Error: imbg not assigned'
    end
    
    if strcmp(fname_conv,'JSO')
        imbg = imbg';  %micromanager images are the inverse of images collected by JSO image collection scripts.
    end
    
    %Convert imbg to double and median filter
    imbg = double(imbg);
    coarse_smooth = 25;
    %smooths background image using course_smooth parameter.  Boundary
    %conditions are symmetric because default 0 bc's causes strange artifacts
    %on the edges.  For these background images symmetric BCs is a good
    %assumption
    imbg = medfilt2(imbg,[coarse_smooth,coarse_smooth],'symmetric');
end
%all locations had 4 sites
%A11 KL 39y1 SDC to SDC+0.5M Sorb
%A12 KL 39y1 SC 2% to no gluc + 0.111M sorb
%B11 KL 39y2 SDC to SDC+0.5M Sorb
%B12 KL 39y2 SC 2% to no gluc + 0.111M sorb
%C11 KL 29y  SDC to SDC+0.5M Sorb
%C12 KL 29y  SC 2% to no gluc + 0.111M sorb
%D11 SC 11-38 SDC to SDC+0.5M Sorb
%D12 SC 11-38 SC 2% to no gluc + 0.111M sorb

legend_vec = {'11-38 0.5M sorb','11-38 No gluc, 0.111M sorb'} %'39y 1 0.5M sorb','39y 1 No gluc, 0.111M sorb','39y 2 0.5M sorb','39y 2 No gluc, 0.111M sorb', '29y 1 0.5M sorb','29y 1 No gluc, 0.111M sorb'}

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc

wellvecSP.SC = {'D11','D12'};
wellvecSP.KL = {'A11','A12','B11','B12'} %RFP: {'A11','A12','B11','B12', 'C11','C12'}; %YFP: {'A11','A12','B11','B12'}
Nsites = 4;

for sp = 1:length(species_cell);
    wellvec = wellvecSP.(species_cell{sp});
    for jj = 1:length(wellvec);
        for kk = 1:Nsites;
            posvecSP.(species_cell{sp}){jj,kk} = [wellvec{jj},'-Site_',num2str(kk-1)];
        end
    end
end

%Remove various sites from list 
%G1 site 2 has no cells
%posvecSP.SC{7,3} = 'NA'
%difficult combinatorics for E1 site 0 - some big cells right next to each
%other but nothing obvious
%posvecSP.SC{5,1} = 'NA'
%posvecSP.SC{5,3} = 'NA'



% posvec.SC = {'A7_site1','p2','p3';
% 'p6','p7','p8';
% 'p11','p12','p13';
% 'p16','p17','p18';
% 'p21','p22','p23';
% 'p26','p27','p28'};
% 
% posvec.KL = {'p4','p5';
% 'p9','p10';
% 'p14','p15';
% 'p19','p20';
% 'p24','p25';
% 'p29','p30'};

%Obtain and store data for each dose (note: only need to do this once) 

get_data = 0;

all_tracks_vec = [];
all_times_vec = [];
posvec = posvecSP.(species);
if get_data == 1 
    for jj = 1:length(wellvecSP.(species))
       for ph = 1:length(phases);
            phase = phases{ph};
            imdir = imdirPhase.(phase);
            time_calc = time_calc_phase.(phase);
            pos_fnames = posvec(jj,:);   %{'p16','p17','p18'} %, 'p14','p15'}
            %remove bad positions (NAs)
            pos_fnames = pos_fnames(~strcmp(pos_fnames,'NA'));
            %this function should take a list of positions as an input
            [tracks,times] = KL_vs_SC_analysis_2color(ipdir,storeim,fname_conv,op_amp,std_threshSP,species, imdir, maxdisp_1x,pos_fnames,channels,channels_to_image,time_calc,imbg);
            all_tracks.(phase) = tracks;
            all_times.(phase) = times + shift_timing(ph);
       end
       all_tracks_vec{jj} = all_tracks;
       all_times_vec{jj} = all_times;
    end
    %store data
    save([base_dir,fname_saveSP.(species)],'all_times_vec','all_tracks_vec','posvec')
else
    %retreive data
    load([base_dir,fname_saveSP.(species)],'all_times_vec','all_tracks_vec','posvec')   
end


%the site I want to image is B11 which is the first site in the stored
%data. A11 has a bad first image.

sites = {'B11','B12'};
site_list = {'A11','A12','B11','B12'};

%normalize by transforming all points linearly from the range spanning the
%minimum to the maximum value of the mean for all experiments to 0-1.
norm_val_max.(channels_to_image{1}) = 1;
norm_val_max.(channels_to_image{2}) = 1;
norm_val_min.(channels_to_image{1}) = 100;
norm_val_min.(channels_to_image{2}) = 100;

for si = 1:length(sites);
    site = sites{si};
    site_ind = find(strcmp(site_list,site));
    all_tracks = all_tracks_vec{site_ind};
    all_times = all_times_vec{site_ind};

    for ph = 1:length(phases)
        figure(ph)
        clf
        phase = phases{ph};
        Ntracks = length(all_tracks.(phase));
        Ntimes = length(all_times.(phase));
        times_ind = 1:Ntimes;
        for ch = 1:length(channels_to_image)
            channel = channels_to_image{ch};
            %Build data storage matrix
            sing_cell_tracks_nfmat = zeros(Ntracks,Ntimes);
             %Go through each track and place it in the correct row at the
             %appropriate timepoint.
            clear tracks
            for jj = 1:Ntracks;
                  tracks_row_nf = [all_tracks.(phase)(jj).nf.(channel)];
                  tracks_row_nmi = [all_tracks.(phase)(jj).nmi.(channel)];
                  tracks_row_times = [all_tracks.(phase)(jj).times];
                  sing_cell_tracks_nfmat(jj,[all_tracks.(phase)(jj).times])= tracks_row_nf;
                  tracks(jj).nf = tracks_row_nf;
                  tracks(jj).nmi = tracks_row_nmi;
                  tracks(jj).times = tracks_row_times;
            end
            
            [mean_nf, std_nf] = nf_calcs(times_ind,tracks);
            sing_cell_tracks.(site).(phase).nf_mat.(channel) = sing_cell_tracks_nfmat;
            sing_cell_tracks.(site).(phase).nf_mean.(channel) = mean_nf;
            norm_val_max.(channel) = max(norm_val_max.(channel),max(mean_nf));
            norm_val_min.(channel) = min(norm_val_min.(channel),min(mean_nf));
            sing_cell_tracks.(site).(phase).nf_std.(channel) = std_nf;
            [mean_nmi, std_nmi] = nmi_calcs(times_ind,tracks);
            sing_cell_tracks.(site).(phase).nmi_mean.(channel) = mean_nmi;
            sing_cell_tracks.(site).(phase).nmi_std.(channel) = std_nmi;            
            sing_cell_tracks.(site).(phase).nmi_std.(channel) = std_nmi;
            sing_cell_tracks.(site).(phase).times = all_times.(phase);
        end
        
    end
    
    
end

%Normalize by norm_val
%plot
phase = 'Post'
site = 'B12'
%condition = '2% Glu -> 0.5M Sorbitol';
condition = '2% Glu -> no gluc, 0.11M Sorb';

figure(1)
channel = channels_to_image{1};
sing_cell_mat_1 = sing_cell_tracks.(site).(phase).nf_mat.(channel);
sing_cell_mat_1_norm = (sing_cell_mat_1-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));

%sort rows: by max value in the first channel, first 4 timepoints.
[sing_cell_mat_1_norm_sorted,sort_ind] = sortrows(sing_cell_mat_1_norm,[-1,-2,-3,-4]);
zero_shift = (0-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
sing_cell_mat_1_norm_sorted(sing_cell_mat_1_norm_sorted==zero_shift) = nan;

[nr,nc] = size(sing_cell_mat_1_norm_sorted);
pcolor([sing_cell_mat_1_norm_sorted nan(nr,1); nan(1,nc+1)]);
colormap(cool)
shading flat;
set(gca, 'ydir', 'reverse');
colorbar
title(['KL: KL.MSN2 nuclear localization. ', condition, '.'])
ylabel('Cells')
xlabel('Time')
set(gca,'XTickLabel',sprintf('%0.2f|',[sing_cell_tracks.(site).(phase).times]))
caxis([-0.3,3.1])

figure(2)
channel = channels_to_image{2};
sing_cell_mat_2 = sing_cell_tracks.(site).(phase).nf_mat.(channel);
sing_cell_mat_2_norm = (sing_cell_mat_2-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
%sort using index from previous channel
sing_cell_mat_2_norm_sorted = sing_cell_mat_2_norm(sort_ind,:);
zero_shift = (0-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
sing_cell_mat_2_norm_sorted(sing_cell_mat_2_norm_sorted==zero_shift) = nan;

[nr,nc] = size(sing_cell_mat_2_norm_sorted);
pcolor([sing_cell_mat_2_norm_sorted nan(nr,1); nan(1,nc+1)]);
colormap(cool)
shading flat;
set(gca, 'ydir', 'reverse');
colorbar
title(['KL: SC.MSN2 nuclear localization. ', condition, '.'])
ylabel('Cells')
xlabel('Time')
set(gca,'XTickLabel',sprintf('%0.2f|',[sing_cell_tracks.(site).(phase).times]))
caxis([-0.3,3.1])

figure(3)
clf 
hold on
cmap = cool(length([sing_cell_tracks.(site).(phase).times]));
[Ncells,Ntimes] = size(sing_cell_mat_1);
for jj = 1:Ntimes;
    cell_vec_1 = sing_cell_mat_1(:,jj);
    cell_vec_1 = cell_vec_1(sing_cell_mat_1(:,jj)>0);
    cell_vec_2 = sing_cell_mat_2(:,jj);
    cell_vec_2 = cell_vec_2(sing_cell_mat_1(:,jj)>0); %used same filter intentionally - although should be the same
    scatter(cell_vec_1,cell_vec_2,5,cmap(jj,:),'filled')
    corr_vec(jj) = corr(cell_vec_1,cell_vec_2);
end
xlabel('KL.MSN2(RFP)')
ylabel('SC.MSN2(YFP)')
title(['KL: ',condition])
axis([1,4,1,4.5])
%scatterplot of data altering color

figure(4)
plot([sing_cell_tracks.(site).(phase).times],corr_vec,'LineWidth',3)
title(['KL: ',condition,'. Correlation KL.MSN2 n.loc. to SC.MSN2 n.loc'])
xlabel('Time')
ylabel('Corr')
axis([10,55,0,1])
return
   
%set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
%cmap = jet(length(legend_vec));
cmap = [1,0,0;
0,0,0];



figure(1)
clf 
hold on
for jj = 1:length(legend_vec)
    all_tracks = all_tracks_vec{jj};
    all_times = all_times_vec{jj};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5);
        set(p,'Parent',plt_grp(jj))
    end
end



hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Condition')
title('MSN2 Nuclear Localization after media change')


figure(2)
clf 
hold on
for jj = 1:length(legend_vec)
    all_tracks = all_tracks_vec{jj};
    all_times = all_times_vec{jj};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0,'nmi','linewidth',1.5);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Condition')
title('MSN2 median intensity')


figure(3)
clf 
hold on
for jj = 1:length(legend_vec)
    all_tracks = all_tracks_vec{jj};
    all_times = all_times_vec{jj};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_ncells(times,tracks,color_val);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Condition')
title('Number of cells identified')



figure(4)
clf 
hold on
legend_vec_RFP = {'(KL.MSN2) 39y 1 No gluc, 0.111M sorb';
    '(KL.MSN2) 39y 1 0.5M sorb';
    %'(KL.MSN2) 39y 2 0.5M sorb';
    %'(KL.MSN2) 39y 2 No gluc, 0.111M sorb';
    '(KL.MSN2) 29y 1 No gluc, 0.111M sorb';
    '(KL.MSN2) 29y 1 0.5M sorb'}

fname_save = '20140703_processed_data_KL_RFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')   
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;

% cmap = [0,0,1; 
% 0,1,0.5; 
% 0,.3333,.8333;
% 0,0.6667,0.6667;
% 1,0,0;
% 0,0,0]

cmap = [0,0,1; 
 .5,0.25,0; 
 %0,.3333,.8333;
 %0,0.6667,0.6667;
 0,0,0;
 1,0,0];

perm = [2,1,6,5];

for jj = 1:length(legend_vec_RFP)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5);
        set(p,'Parent',plt_grp(jj))
    end
end


% legend_vec_YFP = {'(SC.MSN2) 39y 1 0.5M sorb';
%     '(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
%     '(SC.MSN2) 39y 2 0.5M sorb';
%     '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'};

legend_vec_YFP = {'(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
    '(SC.MSN2) 39y 1 0.5M sorb';
%    '(SC.MSN2) 39y 2 0.5M sorb';
%    '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'
     };

cmap = [0,0,1; 
 .5,0.25,0; 
 %0,.3333,.8333;
 %0,0.6667,0.6667;
 %0,0,0;
 %1,0,0
 ];

perm = [2,1];
 
fname_save = '20140703_processed_data_KL_YFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
%all_times_vec_YFP = all_times_vec;
%all_tracks_vec_YFP = all_tracks_vec;
%posvec_YFP = posvec;

for jj = 1:length(legend_vec_YFP)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    kk = jj+length(legend_vec_RFP); %step up plot group
    plt_grp(kk) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5,'LineStyle',':');
        set(p,'Parent',plt_grp(kk))
    end
end

%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 and SC.MSN2 response to glucose dropout and sorbitol in KL')
xlabel('time')
ylabel('Nuclear Localization')



%For a particular condition
site = 1;   %this is the index in the saved data for A11,KL 39y1 SDC to SDC+0.5M Sorb 
%For a particular phase
phase = 'Post'

%load data for each color
fname_save = '20140703_processed_data_KL_RFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
all_times_vec_RFP = all_times_vec;
all_tracks_vec_RFP = all_tracks_vec;

all_times_RFP =  all_times_vec_RFP{site}.(phase);
all_tracks_RFP = all_tracks_vec_RFP{site}.(phase);

fname_save = '20140703_processed_data_KL_YFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
all_times_vec_YFP = all_times_vec;
all_tracks_vec_YFP = all_tracks_vec;

all_times_YFP = all_times_vec_YFP{site}.(phase);
all_tracks_YFP = all_tracks_vec_YFP{site}.(phase);

if length(all_times_YFP) ~= length(all_times_RFP)
    'Error - different number of time points in each channel'
end

times_ind = 1:length(all_times_vec_YFP);

%Find mean 
[nf_mean_RFP, nf_std_RFP] = nf_calcs(times_ind,tracks);
[nf_mean_RFP, nf_std_RFP] = nf_calcs(times_ind,tracks);




profile off
end



