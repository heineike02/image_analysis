function all_tracks = BMH_20140716_analysis_37_38_GD()

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

base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20140716_37_38_GD\'
imdirPhase.Pre = [base_dir,'Pre\']
imdirPhase.Post = [base_dir,'Post\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%input information from experiment here
species = 'KL' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species
species_cell = {'SC','KL'}

channels = {'BF','RFP','YFP'}
channel_to_image = 'YFP'

fname_saveSP.KL = ['20140716_processed_data_KL_', channel_to_image, '.mat'];
fname_saveSP.SC = ['20140716_processed_data_SC_', channel_to_image, '.mat'];

phases =  {'Pre','Post'} %,'Post'} 
shift_timing = [0,8.66666]    
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
%Extract times from metadata in micromanager images by first generating
%metadata_parsed.txt files with python
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
%A8 SC 37_1 GD 
%B8 SC 37_2 GD 
%C8 SC 38_1 GD 
%D8 SC 38_2 GD 
%E8 SC 11-38 GD 
%F8 KL 39 GD  
%G8 SC 37_1 SDC to SDC+0.5M Sorb

%Note: First image after perturbation is at slightly different location
%than second image - shifted about 1/3 of the size of the image.  Can trust
%summary values but not cell tracks.  Should leave out of cell tracking

%SC 
legend_vec =   {'KL 39 GD'}; %{'37_1 GD','37_2 GD','38_1 GD','38_2 GD','11-38 GD','37_1 0.5M Sorb'}; % {'37_1 GD','37_2 GD','37_1 0.5M Sorb'}; %{'KL 39 GD'};
%{'37_1 GD','37_2 GD','38_1 GD','38_2 GD','11-38 GD','37_1 0.5M Sorb'};

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc

wellvecSP.SC = {'A8','B8','C8','D8','E8','G8'};
wellvecSP.KL = {'F8'};
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
    for jj = 1:length(legend_vec)
       for ph = 1:length(phases);
            phase = phases{ph};
            imdir = imdirPhase.(phase);
            time_calc = time_calc_phase.(phase);
            pos_fnames = posvec(jj,:);   %{'p16','p17','p18'} %, 'p14','p15'}
            %remove bad positions (NAs)
            pos_fnames = pos_fnames(~strcmp(pos_fnames,'NA'));
            %this function should take a list of positions as an input
            [tracks,times] = KL_vs_SC_analysis(ipdir,storeim,fname_conv,op_amp,std_threshSP,species, imdir, maxdisp_1x,pos_fnames,channels,channel_to_image,time_calc,imbg);
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

%set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
cmap = jet(length(legend_vec));
%cmap = [1,0,0;
%0,0,0];


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


%return 


figure(4)
clf 
hold on
legend_vec_RFP = {'(SC.MSN2) 37y_1 No gluc, 0.111M sorb';
    '(SC.MSN2) 37y_1 0.5M sorb';
%    '(SC.MSN2) 37y_2 No gluc, 0.111M sorb';
%    '(SC.MSN2) 11-38 No gluc, 0.111M sorb';
    '(SC.MSN2) 38y_1 No gluc, 0.111M sorb';
%    '(SC.MSN2) 38y_2 No gluc, 0.111M sorb'
    }

fname_save = '20140716_processed_data_SC_RFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')   
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;

cmap = [0,0,1; 
.5,0.25,0; 
%0,0.5,1;
%0.5,.5,.5;
0,0,0;
%0.25,0.25,0.25
];

%perm = [1,6,2,5,3,4];
perm = [1,6,3];

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


legend_vec_YFP = {'(KL.MSN2) 37y_1 No gluc, 0.111M sorb';
    '(KL.MSN2) 37y_1 0.5M sorb';
%    '(KL.MSN2) 37y_2 No gluc, 0.111M sorb'
    }

cmap = [0,0,1; 
.5,0.25,0; 
%0,0.5,1
];

%perm = [1,3,2];
perm = [1,3];

fname_saveSP.KL = '20140716_processed_data_SC_YFP.mat';
load([base_dir,fname_saveSP.(species)],'all_times_vec','all_tracks_vec','posvec') 
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
title('KL.MSN2 and SC.MSN2 response to glucose dropout and sorbitol in SC')
xlabel('time')
ylabel('Nuclear Localization')


figure(5)
clf 
hold on
legend_vec_RFP = {'(KL.MSN2) 29y No gluc, 0.111M sorb'};

fname_save = '20140716_processed_data_KL_RFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')   
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;

cmap = [0,0,1];
perm = [1];

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


legend_vec_YFP = {'(SC.MSN2) 29y No gluc, 0.111M sorb'};

cmap = [0,0,1];
perm = [1];

fname_saveSP.KL = '20140716_processed_data_KL_YFP.mat';
load([base_dir,fname_saveSP.(species)],'all_times_vec','all_tracks_vec','posvec') 
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
title('KL.MSN2 and SC.MSN2 response to glucose dropout and sorbitol in SC')
xlabel('time')
ylabel('Nuclear Localization')

profile off
end



