%function [] = Figure_Msn2_Nuc_loc()
%% Set Image Processing directory

%profile off
ipdir = 'C:\Users\BMH_work\github\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)




%% K.Lac NMPP1 for BMH 55, 75, 77 and 78  
data_dir = 'C:\Users\BMH_work\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\'
base_dir = strcat(data_dir, '20151113_55_75_78_prepost_nmpp1\')

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    
x_shift =  8.2  


% Sorb
% RFP Channel is KL.MSN2
%
% Row #	Strain	Media 		YFP channel   AS mutants	Crispr plas
% A1 	BMH55y	4uM 1NMPP1	NA              NA          -
% B1	BMH75y	"       	SC.MSN2         TPK2        -
% C1	BMH77y	"       	"               TPK2/3      +
% D1	BMH78y	"       	"               "           -
% E1	BMH55y	SDC       	NA              NA          -
% F1	BMH75y	"          	SC.MSN2         TPK2        -
% G1	BMH77y	"       	"               TPK2/3      +
% H1	BMH78y	"       	"               TPK2/3      -


fname_save = 'processed_data_75_78.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec_75_78 = all_tracks_vec;
all_times_vec_75_78 = all_times_vec;

fname_save = 'processed_data_55.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec_55 = all_tracks_vec;
all_times_vec_55 = all_times_vec;

all_tracks_vec = [all_tracks_vec_75_78, all_tracks_vec_55];
all_times_vec = [all_times_vec_75_78, all_times_vec_55];

%% RFP channel 

channel = 'RFP'
legend_vec_RFP = {'TPK2-AS, 1-NM-PP1', 'TPK2/3-AS +Cas9, 1-NM-PP1','PKA-AS, 1-NM-PP1','TPK2-AS SDC','TPK2/3-AS +Cas9, SDC','PKA-AS, No Drug','WT, 1-NM-PP1', 'WT, No Drug'}
%B, C, D, F, G, H, A, E
%4uM 1-NM-PP1
%SDC
cmap_RFP = [       0,0,0.0;  %Green
     0,0,0;  %black
     32, 122, 170;  %dark blue
    0,0,0;  %black
    0,0,0;  %black
    32, 122, 170;  %dark blue
    0,0,0;  %black
    0,0,0;  %black
          ]/255;  

linestyle_vec_RFP = {'--','-.','-','--','-.','--','-','--'}
marker_vec_RFP = {'None','None','s','None','None','o','s','o'}

%KL.MSN2 for all strains
figure(1)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [3,7,6,8];
N_RFP = length(perm);

legend_vec_RFP_plot = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plt_grp = zeros(length(legend_vec_RFP_plot),1);
for jj = 1:N_RFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_RFP{perm(jj)},'Marker',marker_vec_RFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph})-x_shift;
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
axis([-5,13,1.2,3])
title('K.Lactis')
xlabel('time(min)')
ylabel('KL.Msn2 Nuclear Localization')

% %% YFP Channel
% 
% fname_save = 'processed_data_75_78.mat';
% load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
% 
% channel = 'YFP'
% legend_vec_YFP = {'TPK2-AS 4uM NMPP1', 'TPK2/3-AS +Cas9 4uM NMPP1','TPK2/3-AS 4uM NMPP1','TPK2-AS SDC','TPK2/3-AS +Cas9 SDC','TPK2/3-AS SDC'}
% 
% cmap_YFP = [       0,0,0.5;  %light Blue
%      0,0,0.5;  %light Blue
%      0,0,0.5;  %light Blue
%     0,0,0;  %black
%     0,0,0;  %black
%     0,0,0  %black
%           ]; 
%       
% linestyle_vec_YFP = {'--','-.','-','--','-.','-'}
% marker_vec_YFP = {'None','None','None','None','None','None'}
% 
% figure(3)
% clf
% hold on
% 
% perm = [2,3,1,5,6,4];
% N_YFP = length(perm);
% 
% legend_vec_YFP_plot = legend_vec_YFP(perm);
% cmap = cmap_YFP(perm,:);
% plt_grp = zeros(length(legend_vec_YFP_plot),1);
% for jj = 1:N_YFP
%     all_tracks = all_tracks_vec{perm(jj)};
%     all_times = all_times_vec{perm(jj)};
%     color_val = cmap(jj,:);
%     plt_grp(jj) = hggroup;
%     plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}, 'Marker',marker_vec_YFP{perm(jj)}};
%     for ph = 1: length(phases)
%         tracks = all_tracks.(phases{ph});
%         timevals = all_times.(phases{ph});
%         p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
%         set(p,'Parent',plt_grp(jj))
%     end
% end
% 
% hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
% axis([0,100,1,5])
% title('SC.Msn2-YFP')
% xlabel('time')
% ylabel('Nuclear Localization')
% 
% 
% figure(4)
% clf
% hold on
% 
% perm = [3,1,6,4];
% N_YFP = length(perm);
% 
% legend_vec_YFP_plot = legend_vec_YFP(perm);
% cmap = cmap_YFP(perm,:);
% plt_grp = zeros(length(legend_vec_YFP_plot),1);
% for jj = 1:N_YFP
%     all_tracks = all_tracks_vec{perm(jj)};
%     all_times = all_times_vec{perm(jj)};
%     color_val = cmap(jj,:);
%     plt_grp(jj) = hggroup;
%     plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}, 'Marker',marker_vec_YFP{perm(jj)}};
%     for ph = 1: length(phases)
%         tracks = all_tracks.(phases{ph});
%         timevals = all_times.(phases{ph});
%         p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
%         set(p,'Parent',plt_grp(jj))
%     end
% end
% 
% hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
% axis([0,100,1,5])
% title('SC.Msn2-YFP')
% xlabel('time')
% ylabel('Nuclear Localization')


%% SC NMPP1 Does Response
base_dir = strcat(data_dir, '20160929_SC_NMPP1_Dose_Resp\')   %20160929_SC_NMPP1_Dose_Resp

phases =  {'Exp'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    

x_shift = 7.5  %6.49

% Want to see 
%
% C2 (AS, 4uM)  9
% D2 (WT, 4uM)  10
% A1 (AS, 0uM)  1
% B1 (WT, 0uM)  2

% 
% Well	Strain	Concentration	Results (visual)
% A1	AS	0	No loc
% B1	WT	0	
% C1	AS	10nM	
% D1	AS	25nM	
% E1	AS	100nM	
% F1	WT	100nM	
% A2	AS	250nM	Some Loc
% B2	AS	1uM	
% C2	AS	4uM	
% D2	WT 	4uM	
% E2	AS	10uM	
% F2	WT	10uM	
% G2	WT	0 10 DMSO: 490 SDC	


%Timeframe to focus on: From addition of drug (aftersecond timepoint (we
%will call that

fname_save = 'processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec;
all_times_vec;

%% RFP channel 

channel = 'RFP'
legend_vec_RFP = {'PKA-AS, No Drug','WT, No Drug', 'AS 10nM', 'AS 25 nM', 'AS 100nM', 'WT 100nM', 'AS 250nM', 'AS 1uM', 'PKA-AS, 1-NM-PP1', 'WT, 1-NM-PP1', 'AS 10uM', 'WT 10uM','WT DMSO'}
cmap_RFP = [ 32, 122, 170;  %dark blue
0, 0, 0;  %black
67, 42, 27;
156, 97, 62;
253, 158, 100;
253, 158, 100;
254, 160, 102;
255, 161, 103;
32, 122, 170;  %dark blue
0, 0, 0; %black
255, 199, 127;
255, 199, 127;
0, 0, 0
]/255;

linestyle_vec_RFP = {'--','--','-','-','-','--','-','-','-','-','-','--','-.'}
marker_vec_RFP = {'o','o','None','None','None','None','None','None','s','s','None','None','None' }

%All Strains and doses
figure(2)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [9,10,1,2];
N_RFP = length(perm);

legend_vec_RFP_plot = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plt_grp = zeros(length(legend_vec_RFP_plot),1);
for jj = 1:N_RFP
    all_tracks = all_tracks_vec{perm(jj)} ;
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_RFP{perm(jj)},'Marker',marker_vec_RFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph})-x_shift;
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
axis([-5,13,1,9])
title('S.Cerevisiae')
xlabel('time(min)')
ylabel('Msn2 Nuclear Localiztion')  %(\Delta-ZF)-RFP



% %% YFP Channel
% 
% channel = 'YFP'
% legend_vec_YFP = legend_vec_RFP
% cmap_YFP = cmap_RFP
% 
% linestyle_vec_YFP = linestyle_vec_RFP
% marker_vec_YFP = marker_vec_RFP
% 
% %SC.MSN2 for all strains
% figure(2)
% clf
% hold on
% 
% %cmap = jet(length(legend_vec_RFP));
% %perm = 1:length(legend_vec_RFP);
% perm = [2,4,6,1,3,5];
% N_YFP = length(perm);
% 
% legend_vec_YFP_plot = legend_vec_YFP(perm);
% cmap = cmap_YFP(perm,:);
% plt_grp = zeros(length(legend_vec_YFP_plot),1);
% for jj = 1:N_YFP
%     all_tracks = all_tracks_vec{perm(jj)};
%     all_times = all_times_vec{perm(jj)};
%     color_val = cmap(jj,:);
%     plt_grp(jj) = hggroup;
%     plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)},'Marker',marker_vec_YFP{perm(jj)}};
%     for ph = 1: length(phases)
%         tracks = all_tracks.(phases{ph});
%         timevals = all_times.(phases{ph});
%         p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
%         set(p,'Parent',plt_grp(jj))
%     end
% end
% 
% hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
% axis([0,130,1,3])
% title('SC.MSN2-YFP in KL Cells')
% xlabel('time')
% ylabel('Nuclear Localization')
% 
% 
% 
% fname_save = 'processed_data_SC.mat';
% load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
% 
% all_tracks_vec;
% all_times_vec;


% %% SC - from all my pictures: 
% 
% 
% 
% 
% % Filename convention
% % 'JSO','Micromanager', or 'HCS_Nikon'
% % For example of Micromanager filename convention see BMH_20140127 analysis
% % files
% fname_conv = 'Micromanager'
% 
% %imdirPhase.Pre = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_pre_1\'
% %imdirPhase.Post = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_post_1\'
% 
% base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20140918_AS_4uM_p2uM\'
% imdirPhase.Pre = [base_dir,'Pre\']
% imdirPhase.Post = [base_dir,'Post\']
% %imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']
% 
% %input information from experiment here
% species = 'SC' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species
% species_cell = {'SC','KL'} 
% 
% channels = {'BF','RFP','YFP'}
% %too much motion to link cells
% channel_to_image = {'RFP','YFP'}  %if this is just one channel just list it as a text variable i.e. 'RFP'. 
% 
% fname_saveSP.SC = '20140918_processed_data_SC_2Color.mat';
% %fname_saveSP.SC = '20140703_processed_data_SC.mat'; 
% 
% phases =  {'Pre','Post'} %,'Post'} 
% shift_timing = [0,12]    
% %These are the absolute times at which each phase starts.
% %timestep method
% %time_calc:  Tells the program how to calculate each time value.  If the
% %input is a number then that is interpreted as a dt between images
% %If the input is a filename, then that is interpreted as metadata that gives
% %exact time values.  
% 
% %time_calc_phase.Pre = 5
% %time_calc_phase.Post = 4
% %time_calc_phase.Pre = [imdirPhase.Pre, 'acqdat.txt']
% %time_calc_phase.Post = [imdirPhase.Post, 'acqdat.txt']
% %extract actual times from metadata in micromanager images
% %addpath('C:/Users/Ben/Documents/GitHub/image_analysis/jsonlab');
% %time_calc_phase.Pre  = 'metadata.txt'
% %time_calc_phase.Post = 'metadata.txt'
% generate_metadata_parsed = 0;
% metadata_conv_fname = 'C:\Users\Ben\Documents\GitHub\image_analysis\times_from_umanager_metadata.py';
% if generate_metadata_parsed ==1;
%     for ph = [1:length(phases)]
%           phase = phases{ph};
%           imdir = imdirPhase.(phase);
%           system(['python ', metadata_conv_fname, ' ', imdir])
%     end
% end
% time_calc_phase.Pre =  'metadata_parsed.txt';
% time_calc_phase.Post = 'metadata_parsed.txt';
% 
% 
% %std thresh for calling a peak to find a cell
% %0.16 seemed to work best for k.lactis on 16JUL images
% std_threshSP.KL = 0.16;
% std_threshSP.SC = 0.2;
% 
% %Tracking parameters
% %estimated maximum number of tracks
% maxdisp_1x = 4;
% 
% % optical amplification 
% % 1x or 1.5x
% op_amp =  '1.5x'
% storeim = 1;
% 
% bgimg = 1; 
% %1 if you have a background image, 0 if you don't. 
% %Gets background image depending on channel to image
% %Collect background images using micromanager 
% if bgimg == 0
%     imbg = 1; % default
% else
%     %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
%     if iscell(channel_to_image)    
%         ch2i = channel_to_image;
%     else
%         ch2i = {channel_to_image};
%     end
%    
%     for jj = 1:length(ch2i)
%         ch2i_txt = ch2i{jj}
%         if strcmp(ch2i_txt,'RFP')
%              imbg_jj = imread([base_dir,'BG\img_000000000_RFP_001.tif']);
%         elseif strcmp(ch2i_txt,'YFP')
%              imbg_jj = imread([base_dir,'BG\img_000000000_YFP_001.tif']);
%             %imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\YFP_BG\img_000000000_Default_000.tiff');
%         else
%              'Error: imbg not assigned'
%         end
%     
%         if strcmp(fname_conv,'JSO')
%             imbg_jj = imbg_jj';  %micromanager images are the inverse of images collected by JSO image collection scripts.
%         end
%     
%         %Convert imbg to double and median filter
%         imbg_jj = double(imbg_jj);
%         coarse_smooth = 25;
%         %smooths background image using course_smooth parameter.  Boundary
%         %conditions are symmetric because default 0 bc's causes strange artifacts
%         %on the edges.  For these background images symmetric BCs is a good
%         %assumption
%         imbg_jj = medfilt2(imbg_jj,[coarse_smooth,coarse_smooth],'symmetric');
%         imbg.(ch2i_txt) = imbg_jj;
% end
% %all locations had 4 sites
% %4uM 1-NM-PP1 
% %A1 37, 
% %A2 42-1, 
% %A3 42-2, 
% %0.2uM: 
% %A4 42-1, 
% %A5 42-2 
% %0 
% %A6 37 
% %A7 42-1, 
% %A8 42-2
% 
% legend_vec = {'37 4uM','42_1 4uM', '42_2 4uM','42_1 0.2 uM1' , '42_2 0.2 uM', '37 cont','42_1 cont', '42_2 cont'} %'39y 1 0.5M sorb','39y 1 No gluc, 0.111M sorb','39y 2 0.5M sorb','39y 2 No gluc, 0.111M sorb', '29y 1 0.5M sorb','29y 1 No gluc, 0.111M sorb'}
% 
% %Micromanager positions: Pos0, Pos1, etc.  
% %JSO positions p1,p2,etc
% 
% wellvecSP.SC = {'A1','B1','C1','D1', 'E1', 'F1', 'G1', 'H1'};  %RFP_only 'C8','D8','E8'
% wellvecSP.KL = {'F8'};
% Nsites = 4;
% 
% for sp = 1:length(species_cell);
%     wellvec = wellvecSP.(species_cell{sp});
%     for jj = 1:length(wellvec);
%         for kk = 1:Nsites;
%             posvecSP.(species_cell{sp}){jj,kk} = [wellvec{jj},'-Site_',num2str(kk-1)];
%         end
%     end
% end
% 
% %Remove various sites from list 
% %C1 site 3 starts off with no cells - Remove
% %B1 site 3 starts off with no cells - Remove
% %A1 Site 3 starts off with one cell
% posvecSP.SC{2,4} = 'NA'
% posvecSP.SC{3,4} = 'NA'
% %difficult combinatorics for E1 site 0 - some big cells right next to each
% %other but nothing obvious
% %posvecSP.SC{5,1} = 'NA'
% %posvecSP.SC{5,3} = 'NA'
% 
% 
% 
% % posvec.SC = {'A7_site1','p2','p3';
% % 'p6','p7','p8';
% % 'p11','p12','p13';
% % 'p16','p17','p18';
% % 'p21','p22','p23';
% % 'p26','p27','p28'};
% % 
% % posvec.KL = {'p4','p5';
% % 'p9','p10';
% % 'p14','p15';
% % 'p19','p20';
% % 'p24','p25';
% % 'p29','p30'};
% 
% %Obtain and store data for each dose (note: only need to do this once) 
% 
% get_data = 0;
% 
% all_tracks_vec = [];
% all_times_vec = [];
% posvec = posvecSP.(species);
% if get_data == 1 
%     for jj = 1:length(wellvecSP.(species))
%        for ph = 1:length(phases);
%             phase = phases{ph};
%             imdir = imdirPhase.(phase);
%             time_calc = time_calc_phase.(phase);
%             pos_fnames = posvec(jj,:);   %{'p16','p17','p18'} %, 'p14','p15'}
%             %remove bad positions (NAs)
%             pos_fnames = pos_fnames(~strcmp(pos_fnames,'NA'));
%             %this function should take a list of positions as an input
%             [tracks,times] = KL_vs_SC_analysis(ipdir,storeim,fname_conv,op_amp,std_threshSP,species, imdir, maxdisp_1x,pos_fnames,channels,channel_to_image,time_calc,imbg);
%             all_tracks.(phase) = tracks;
%             all_times.(phase) = times + shift_timing(ph);
%        end
%        all_tracks_vec{jj} = all_tracks;
%        all_times_vec{jj} = all_times;
%     end
%     %store data
%     save([base_dir,fname_saveSP.(species)],'all_times_vec','all_tracks_vec','posvec')
% else
%     %retreive data
%     load([base_dir,fname_saveSP.(species)],'all_times_vec','all_tracks_vec','posvec')   
% end
% 
% 
% %the site I want to image is B11 which is the first site in the stored
% %data. A11 has a bad first image.
% 
% sites ={'A1','B1','C1','D1', 'E1', 'F1', 'G1', 'H1'};
% wellvecSP.SC = {'A1','B1','C1','D1', 'E1', 'F1', 'G1', 'H1'};
% 
% %normalize by transforming all points linearly from the range spanning the
% %minimum to the maximum value of the mean for all experiments to 0-1.
% norm_val_max.(channels_to_image{1}) = 1;
% norm_val_max.(channels_to_image{2}) = 1;
% norm_val_min.(channels_to_image{1}) = 100;
% norm_val_min.(channels_to_image{2}) = 100;
% 
% for si = 1:length(sites);
%     site = sites{si};
%     site_ind = find(strcmp(site_list,site));
%     all_tracks = all_tracks_vec{site_ind};
%     all_times = all_times_vec{site_ind};
% 
%     for ph = 1:length(phases)
%         figure(ph)
%         clf
%         phase = phases{ph};
%         Ntracks = length(all_tracks.(phase));
%         Ntimes = length(all_times.(phase));
%         times_ind = 1:Ntimes;
%         for ch = 1:length(channels_to_image)
%             channel = channels_to_image{ch};
%             %Build data storage matrix
%             sing_cell_tracks_nfmat = zeros(Ntracks,Ntimes);
%              %Go through each track and place it in the correct row at the
%              %appropriate timepoint.
%             clear tracks
%             for jj = 1:Ntracks;
%                   tracks_row_nf = [all_tracks.(phase)(jj).nf.(channel)];
%                   tracks_row_nmi = [all_tracks.(phase)(jj).nmi.(channel)];
%                   tracks_row_times = [all_tracks.(phase)(jj).times];
%                   sing_cell_tracks_nfmat(jj,[all_tracks.(phase)(jj).times])= tracks_row_nf;
%                   tracks(jj).nf = tracks_row_nf;
%                   tracks(jj).nmi = tracks_row_nmi;
%                   tracks(jj).times = tracks_row_times;
%             end
%             
%             [mean_nf, std_nf] = nf_calcs(times_ind,tracks);
%             sing_cell_tracks.(site).(phase).nf_mat.(channel) = sing_cell_tracks_nfmat;
%             sing_cell_tracks.(site).(phase).nf_mean.(channel) = mean_nf;
%             norm_val_max.(channel) = max(norm_val_max.(channel),max(mean_nf));
%             norm_val_min.(channel) = min(norm_val_min.(channel),min(mean_nf));
%             sing_cell_tracks.(site).(phase).nf_std.(channel) = std_nf;
%             [mean_nmi, std_nmi] = nmi_calcs(times_ind,tracks);
%             sing_cell_tracks.(site).(phase).nmi_mean.(channel) = mean_nmi;
%             sing_cell_tracks.(site).(phase).nmi_std.(channel) = std_nmi;            
%             sing_cell_tracks.(site).(phase).nmi_std.(channel) = std_nmi;
%             sing_cell_tracks.(site).(phase).times = all_times.(phase);
%         end
%         
%     end
%     
%     
% end
% 
% 
% %{
% %Normalize by norm_val
% %plot
% phase = 'Post'
% site = 'A8'
% %condition = '2% Glu -> 0.5M Sorbitol';
% condition = '2% Glu -> no gluc, 0.11M Sorb';
% 
% figure(1)
% channel = channels_to_image{1};
% sing_cell_mat_1 = sing_cell_tracks.(site).(phase).nf_mat.(channel);
% sing_cell_mat_1_norm = (sing_cell_mat_1-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
% 
% %sort rows: by max value in the first channel, first 4 timepoints.
% [sing_cell_mat_1_norm_sorted,sort_ind] = sortrows(sing_cell_mat_1_norm,[-1,-2,-3,-4]);
% zero_shift = (0-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
% sing_cell_mat_1_norm_sorted(sing_cell_mat_1_norm_sorted==zero_shift) = nan;
% 
% [nr,nc] = size(sing_cell_mat_1_norm_sorted);
% pcolor([sing_cell_mat_1_norm_sorted nan(nr,1); nan(1,nc+1)]);
% colormap(cool)
% shading flat;
% set(gca, 'ydir', 'reverse');
% colorbar
% title(['SC: SC.MSN2 nuclear localization. ', condition, '.'])
% ylabel('Cells')
% xlabel('Time')
% set(gca,'XTick',1:length([sing_cell_tracks.(site).(phase).times])+0.5)
% set(gca,'XTickLabel',sprintf('%0.0f|',[sing_cell_tracks.(site).(phase).times]))
% caxis([-0.4,4])
% 
% figure(2)
% channel = channels_to_image{2};
% sing_cell_mat_2 = sing_cell_tracks.(site).(phase).nf_mat.(channel);
% sing_cell_mat_2_norm = (sing_cell_mat_2-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
% %sort using index from previous channel
% sing_cell_mat_2_norm_sorted = sing_cell_mat_2_norm(sort_ind,:);
% zero_shift = (0-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
% sing_cell_mat_2_norm_sorted(sing_cell_mat_2_norm_sorted==zero_shift) = nan;
% 
% [nr,nc] = size(sing_cell_mat_2_norm_sorted);
% pcolor([sing_cell_mat_2_norm_sorted nan(nr,1); nan(1,nc+1)]);
% colormap(cool)
% shading flat;
% set(gca, 'ydir', 'reverse');
% colorbar
% title(['SC: KL.MSN2 nuclear localization. ', condition, '.'])
% ylabel('Cells')
% xlabel('Time')
% set(gca,'XTick',1:length([sing_cell_tracks.(site).(phase).times])+0.5)
% set(gca,'XTickLabel',sprintf('%0.0f|',[sing_cell_tracks.(site).(phase).times]))
% caxis([-0.4,4])
% 
% figure(3)
% clf 
% hold on
% cmap = cool(length([sing_cell_tracks.(site).(phase).times]));
% [Ncells,Ntimes] = size(sing_cell_mat_1);
% for jj = 1:Ntimes;
%     cell_vec_1 = sing_cell_mat_1(:,jj);
%     cell_vec_1 = cell_vec_1(sing_cell_mat_1(:,jj)>0);
%     cell_vec_2 = sing_cell_mat_2(:,jj);
%     cell_vec_2 = cell_vec_2(sing_cell_mat_1(:,jj)>0); %used same filter intentionally - although should be the same
%     scatter(cell_vec_1,cell_vec_2,5,cmap(jj,:),'filled')
%     corr_vec(jj) = corr(cell_vec_1,cell_vec_2);
% end
% xlabel('SC.MSN2(RFP)')
% ylabel('KL.MSN2(YFP)')
% title(['SC: ',condition])
% axis([1,14,1,9])
% %scatterplot of data altering color
% 
% figure(4)
% plot([sing_cell_tracks.(site).(phase).times],corr_vec,'LineWidth',3)
% title(['SC: ',condition,'. Correlation SC.MSN2 n.loc. to KL.MSN2 n.loc'])
% xlabel('Time')
% ylabel('Corr')
% axis([10,70,0,1])
% 
% %set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
% %cmap = jet(length(legend_vec));
% cmap = [1,0,0;
% 0,0,0];
% 
% %}
% 
% figure(1)
% clf 
% hold on
% for jj = 1:length(legend_vec)
%     all_tracks = all_tracks_vec{jj};
%     all_times = all_times_vec{jj};
%     color_val = cmap(jj,:);
%     plt_grp(jj) = hggroup;
%     for ph = 1: length(phases)
%         tracks = all_tracks.(phases{ph});
%         times = all_times.(phases{ph});
%         p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5);
%         set(p,'Parent',plt_grp(jj))
%     end
% end
% 
% 
% 
% hleg = legend(plt_grp,legend_vec) %,'Location','NE');
% htitle = get(hleg,'Title');
% set(htitle,'String','Condition')
% title('MSN2 Nuclear Localization after media change')
% 
% 
% figure(2)
% clf 
% hold on
% for jj = 1:length(legend_vec)
%     all_tracks = all_tracks_vec{jj};
%     all_times = all_times_vec{jj};
%     color_val = cmap(jj,:);
%     plt_grp(jj) = hggroup;
%     for ph = 1: length(phases)
%         tracks = all_tracks.(phases{ph});
%         times = all_times.(phases{ph});
%         p = plot_meanvalues(times,tracks,color_val,0,'nmi','linewidth',1.5);
%         set(p,'Parent',plt_grp(jj))
%     end
% end
% 
% hleg = legend(plt_grp,legend_vec) %,'Location','NE');
% htitle = get(hleg,'Title');
% set(htitle,'String','Condition')
% title('MSN2 median intensity')
% 
% 
% figure(3)
% clf 
% hold on
% for jj = 1:length(legend_vec)
%     all_tracks = all_tracks_vec{jj};
%     all_times = all_times_vec{jj};
%     color_val = cmap(jj,:);
%     plt_grp(jj) = hggroup;
%     for ph = 1: length(phases)
%         tracks = all_tracks.(phases{ph});
%         times = all_times.(phases{ph});
%         p = plot_ncells(times,tracks,color_val);
%         set(p,'Parent',plt_grp(jj))
%     end
% end
% 
% hleg = legend(plt_grp,legend_vec); %,'Location','NE');
% htitle = get(hleg,'Title');
% set(htitle,'String','Condition')
% title('Number of cells identified')
% 
% 
% 
% figure(4)
% clf 
% hold on
% legend_vec_RFP = {'(KL.MSN2) 39y 1 No gluc, 0.111M sorb';
%     '(KL.MSN2) 39y 1 0.5M sorb';
%     %'(KL.MSN2) 39y 2 0.5M sorb';
%     %'(KL.MSN2) 39y 2 No gluc, 0.111M sorb';
%     '(KL.MSN2) 29y 1 No gluc, 0.111M sorb';
%     '(KL.MSN2) 29y 1 0.5M sorb'}
% 
% fname_save = '20140703_processed_data_KL_RFP.mat';
% load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')   
% %all_times_vec_RFP = all_times_vec;
% %all_tracks_vec_RFP = all_tracks_vec;
% %posvec_RFP = posvec;
% 
% % cmap = [0,0,1; 
% % 0,1,0.5; 
% % 0,.3333,.8333;
% % 0,0.6667,0.6667;
% % 1,0,0;
% % 0,0,0]
% 
% cmap = [0,0,1; 
%  .5,0.25,0; 
%  %0,.3333,.8333;
%  %0,0.6667,0.6667;
%  0,0,0;
%  1,0,0];
% 
% perm = [2,1,6,5];
% 
% for jj = 1:length(legend_vec_RFP)
%     all_tracks = all_tracks_vec{perm(jj)};
%     all_times = all_times_vec{perm(jj)};
%     color_val = cmap(jj,:);
%     plt_grp(jj) = hggroup;
%     for ph = 1: length(phases)
%         tracks = all_tracks.(phases{ph});
%         times = all_times.(phases{ph});
%         p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5);
%         set(p,'Parent',plt_grp(jj))
%     end
% end
% 
% 
% % legend_vec_YFP = {'(SC.MSN2) 39y 1 0.5M sorb';
% %     '(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
% %     '(SC.MSN2) 39y 2 0.5M sorb';
% %     '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'};
% 
% legend_vec_YFP = {'(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
%     '(SC.MSN2) 39y 1 0.5M sorb';
% %    '(SC.MSN2) 39y 2 0.5M sorb';
% %    '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'
%      };
% 
% cmap = [0,0,1; 
%  .5,0.25,0; 
%  %0,.3333,.8333;
%  %0,0.6667,0.6667;
%  %0,0,0;
%  %1,0,0
%  ];
% 
% perm = [2,1];
%  
% fname_save = '20140703_processed_data_KL_YFP.mat';
% load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
% %all_times_vec_YFP = all_times_vec;
% %all_tracks_vec_YFP = all_tracks_vec;
% %posvec_YFP = posvec;
% 
% for jj = 1:length(legend_vec_YFP)
%     all_tracks = all_tracks_vec{perm(jj)};
%     all_times = all_times_vec{perm(jj)};
%     color_val = cmap(jj,:);
%     kk = jj+length(legend_vec_RFP); %step up plot group
%     plt_grp(kk) = hggroup;
%     for ph = 1: length(phases)
%         tracks = all_tracks.(phases{ph});
%         times = all_times.(phases{ph});
%         p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5,'LineStyle',':');
%         set(p,'Parent',plt_grp(kk))
%     end
% end
% 
% %combine legend vectors
% legend_vec = [legend_vec_RFP;legend_vec_YFP];
% hleg = legend(plt_grp,legend_vec) %,'Location','NE');
% htitle = get(hleg,'Title');
% %set(htitle,'String','Condition')
% title('KL.MSN2 and SC.MSN2 response to glucose dropout and sorbitol in KL')
% xlabel('time')
% ylabel('Nuclear Localization')
% 
% 
% 
% %For a particular condition
% site = 1;   %this is the index in the saved data for A11,KL 39y1 SDC to SDC+0.5M Sorb 
% %For a particular phase
% phase = 'Post'
% 
% %load data for each color
% fname_save = '20140703_processed_data_KL_RFP.mat';
% load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
% all_times_vec_RFP = all_times_vec;
% all_tracks_vec_RFP = all_tracks_vec;
% 
% all_times_RFP =  all_times_vec_RFP{site}.(phase);
% all_tracks_RFP = all_tracks_vec_RFP{site}.(phase);
% 
% fname_save = '20140703_processed_data_KL_YFP.mat';
% load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
% all_times_vec_YFP = all_times_vec;
% all_tracks_vec_YFP = all_tracks_vec;
% 
% all_times_YFP = all_times_vec_YFP{site}.(phase);
% all_tracks_YFP = all_tracks_vec_YFP{site}.(phase);
% 
% if length(all_times_YFP) ~= length(all_times_RFP)
%     'Error - different number of time points in each channel'
% end
% 
% times_ind = 1:length(all_times_vec_YFP);
% 
% %Find mean 
% [nf_mean_RFP, nf_std_RFP] = nf_calcs(times_ind,tracks);
% [nf_mean_RFP, nf_std_RFP] = nf_calcs(times_ind,tracks);
% 
% 
% 
% 
% profile off
% end
% 
% 
% 
% 
% 
% 
% 
% 
% return
% 
% 



