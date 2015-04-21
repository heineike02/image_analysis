function all_tracks = BMH_20150408_MSN4_GD_Tm()
%

profile off
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

base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150408_MSN4_GD_TM\'
imdirPhase.Pre = [base_dir,'Pre',filesep]
imdirPhase.Post = [base_dir,'Post',filesep]
%imdirPhase.Post_p2 = [base_dir,'Post_p2\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%input information from experiment here
species = 'SC' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species
species_cell = {'SC','KL'} 

channels = {'BF','RFP','YFP'}
%too much motion to link cells
channel_to_image = {'RFP','YFP'}  %if this is just one channel just list it as a text variable i.e. 'RFP'. 

fname_saveSP.SC = '20150408_processed_data_SC.mat';
fname_saveSP.KL = '20141210_processed_data_KL.mat';
%fname_saveSP.SC = '20140703_processed_data_SC.mat'; 

phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0,8]    
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
metadata_conv_fname = [ipdir,'times_from_umanager_metadata.py'];
if generate_metadata_parsed ==1;
    for ph = [1:length(phases)]
          phase = phases{ph};
          imdir = imdirPhase.(phase);
          phase = ['python ', metadata_conv_fname, ' ', imdir]
          system(['python ', metadata_conv_fname, ' ', imdir])
    end
end
time_calc_phase.Pre =  'metadata_parsed.txt';
time_calc_phase.Post = 'metadata_parsed.txt';
%time_calc_phase.Post_p2 = 'metadata_parsed.txt';

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

bgimg = 1; 
%1 if you have a background image, 0 if you don't. 
%Gets background image depending on channel to image
%Collect background images using micromanager 
if bgimg == 0  %set default to 1
    if iscell(channel_to_image)    
        ch2i = channel_to_image;
    else
        ch2i = {channel_to_image};
    end
    
    for jj = 1:length(ch2i)
        %ch2i_txt = ch2i{jj}
        imbg.(ch2i{jj}) = 1;
    end
    
else
    %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
    if iscell(channel_to_image)    
        ch2i = channel_to_image;
    else
        ch2i = {channel_to_image};
    end
   
    for jj = 1:length(ch2i)
        ch2i_txt = ch2i{jj}
        if strcmp(ch2i_txt,'RFP')
             imbg_jj = imread([base_dir,'BG',filesep,'img_000000000_RFP_001.tif']);
        elseif strcmp(ch2i_txt,'YFP')
             imbg_jj = imread([base_dir,'BG',filesep,'img_000000000_YFP_001.tif']);
            %imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\YFP_BG\img_000000000_Default_000.tiff');
        else
             'Error: imbg not assigned'
        end
    
        if strcmp(fname_conv,'JSO')
            imbg_jj = imbg_jj';  %micromanager images are the inverse of images collected by JSO image collection scripts.
        end
    
        %Convert imbg to double and median filter
        imbg_jj = double(imbg_jj);
        coarse_smooth = 25;
        %smooths background image using course_smooth parameter.  Boundary
        %conditions are symmetric because default 0 bc's causes strange artifacts
        %on the edges.  For these background images symmetric BCs is a good
        %assumption
        imbg_jj = medfilt2(imbg_jj,[coarse_smooth,coarse_smooth],'symmetric');
        imbg.(ch2i_txt) = imbg_jj;
    end
end
%all locations had 4 sites
%A3 - 46 1 MSN4 (TPK AS). Gluc Drop
%B3 - 46 2 MSN4 (TPK AS). Gluc Drop
%C3 - 47 1 MSN4 dZF (TPK AS). Gluc Drop
%D3 - 47 2 MSN4 dZF (TPK AS). Gluc Drop
%E3 - 49 MSN4. Gluc Drop
%F3 - 42 MSN2, KL.MSN2. Gluc Drop
%G3 - 47 SDC
%H3 - 47 Tm 2.5 ug/ml



wellvecSP.SC = {'A3','B3','C3','D3','E3','F3','G3','H3'}; 
wellvecSP.KL = {};
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

%bad data
%C3 had large non-yeast object
posvecSP.SC{3,4} = 'NA'
%posvecSP.SC{6,2} = 'NA'
%posvecSP.SC{6,3} = 'NA'


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


%% 
make_plots = 1

%all locations had 4 sites
%A3 - 46 1 MSN4 (TPK AS). Gluc Drop
%B3 - 46 2 MSN4 (TPK AS). Gluc Drop
%C3 - 47 1 MSN4 dZF (TPK AS). Gluc Drop
%D3 - 47 2 MSN4 dZF (TPK AS). Gluc Drop
%E3 - 49 MSN4. Gluc Drop
%F3 - 42 MSN2, KL.MSN2. Gluc Drop
%G3 - 47 SDC
%H3 - 47 Tm 2.5 ug/ml


if make_plots == 1
    
    %Mean NF time traces MSN2 and MSN4 separately
    figure(1)
    clf
    hold on
    channel = 'RFP'
    legend_vec_RFP = {'(46 SC.MSN2 1) SDC -> SD-Gluc',
        '(46 SC.MSN2 2) SDC -> SD-Gluc',
        '(47 SC.MSN2 1) SDC -> SD-Gluc',
        '(47 SC.MSN2 2) SDC -> SD-Gluc',
        '(49 SC.MSN2) SDC -> SD-Gluc',
        '(42 SC.MSN2) SDC -> SD-Gluc',
        '(47 SC.MSN2 1) SDC -> SDC',
        '(47 SC.MSN2 1) SDC -> Tm 2.5 ug/ml'}
    fname_save = '20150408_processed_data_SC.mat';
    load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
    cmap_RFP = [1,0,0;  %Red
        1,0.35,0; %light red
        0,0,1;  %Blue
        0,0.5,1;  %Light blue
        0,1,0; %green
        0,1,0.5;
        0.2,0.2,0.2; %gray
        0,0,0   %black
        ];
    
    %cmap = jet(length(legend_vec_RFP));
    %perm = 1:length(legend_vec_RFP);
    perm = [1,2,3,4,5,6,7];
    N_RFP = length(perm);
    legend_vec = legend_vec_RFP(perm);
    cmap = cmap_RFP(perm,:);
    plot_params = {'linewidth',1.5,'LineStyle','-'};
    plt_grp = zeros(length(legend_vec),1);
    for jj = 1:length(perm)
        all_tracks = all_tracks_vec{perm(jj)};
        all_times = all_times_vec{perm(jj)};
        color_val = cmap(jj,:);
        plt_grp(jj) = hggroup;
        for ph = 1: length(phases)
            tracks = all_tracks.(phases{ph});
            timevals = all_times.(phases{ph});
            p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
            set(p,'Parent',plt_grp(jj))
        end
    end
    
    hleg = legend(plt_grp,legend_vec); %,'Location','NE');
    htitle = get(hleg,'Title');
    %set(htitle,'String','Condition')
    title('SC Cells, RFP Channel')
    xlabel('time')
    ylabel('Nuclear Localization')
    
    
    figure(2) 
    %Plot NF for just TM cells
    clf
    hold on
    
    %cmap = jet(length(legend_vec_RFP));
    %perm = 1:length(legend_vec_RFP);
    perm = [3,4,7,8]  ;
    N_RFP = length(perm);
    legend_vec = legend_vec_RFP(perm);
    cmap = cmap_RFP(perm,:);
    plot_params = {'linewidth',1.5,'LineStyle','-'};
    plt_grp = zeros(length(legend_vec),1);
    for jj = 1:length(perm)
        all_tracks = all_tracks_vec{perm(jj)};
        all_times = all_times_vec{perm(jj)};
        color_val = cmap(jj,:);
        plt_grp(jj) = hggroup;
        for ph = 1: length(phases)
            tracks = all_tracks.(phases{ph});
            timevals = all_times.(phases{ph});
            p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
            set(p,'Parent',plt_grp(jj))
        end
    end
    
    hleg = legend(plt_grp,legend_vec); %,'Location','NE');
    htitle = get(hleg,'Title');
    %set(htitle,'String','Condition')
    title('SC Cells, RFP Channel')
    xlabel('time')
    ylabel('Nuclear Localization')
    
    
    figure(3)
    
    %Plot NMI for just TM cells
    clf
    hold on
    
    %cmap = jet(length(legend_vec_RFP));
    %perm = 1:length(legend_vec_RFP);
    perm = [3,4,7,8]  ;
    N_RFP = length(perm);
    legend_vec = legend_vec_RFP(perm);
    cmap = cmap_RFP(perm,:);
    plot_params = {'linewidth',1.5,'LineStyle','-'};
    plt_grp = zeros(length(legend_vec),1);
    for jj = 1:length(perm)
        all_tracks = all_tracks_vec{perm(jj)};
        all_times = all_times_vec{perm(jj)};
        color_val = cmap(jj,:);
        plt_grp(jj) = hggroup;
        for ph = 1: length(phases)
            tracks = all_tracks.(phases{ph});
            timevals = all_times.(phases{ph});
            p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nmi','plot_params',plot_params);
            set(p,'Parent',plt_grp(jj))
        end
    end
    
    hleg = legend(plt_grp,legend_vec); %,'Location','NE');
    htitle = get(hleg,'Title');
    %set(htitle,'String','Condition')
    title('SC Cells, RFP Channel')
    xlabel('time')
    ylabel('Median Intensity')
    
    
    
    figure(4)
    clf
    hold on
    channel = 'YFP'
    legend_vec_YFP = {'(46 SC.MSN4 1) SDC -> SD-Gluc',
        '(46 SC.MSN4 2) SDC -> SD-Gluc',
        '(47 SC.MSN4(dzf) 1) SDC -> SD-Gluc',
        '(47 SC.MSN4(dzf) 2) SDC -> SD-Gluc',
        '(49 SC.MSN4(dzf)) SDC -> SD-Gluc',
        '(42 KL.MSN2) SDC -> SD-Gluc',
        '(47 SC.MSN4 1) SDC -> SDC',
        '(47 SC.MSN4 1) SDC -> Tm 2.5 ug/ml'}
    fname_save = '20150408_processed_data_SC.mat';
    load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
    cmap_YFP = cmap_RFP;
    
    %cmap = jet(length(legend_vec_RFP));
    %perm = 1:length(legend_vec_RFP);
    perm = [1,2,3,4,5,6,7] ;
    legend_vec = legend_vec_YFP(perm);
    cmap = cmap_YFP(perm,:);
    plot_params = {'linewidth',1.5,'LineStyle','-'};
    plt_grp = zeros(length(legend_vec),1);
    for jj = 1:length(perm)
        all_tracks = all_tracks_vec{perm(jj)};
        all_times = all_times_vec{perm(jj)};
        color_val = cmap(jj,:);
        plt_grp(jj) = hggroup;
        for ph = 1: length(phases)
            tracks = all_tracks.(phases{ph});
            timevals = all_times.(phases{ph});
            p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
            set(p,'Parent',plt_grp(jj))
        end
    end
    
    
    hleg = legend(plt_grp,legend_vec) %,'Location','NE');
    htitle = get(hleg,'Title');
    %set(htitle,'String','Condition')
    title('SC Cells, YFP Channel')
    xlabel('time')
    ylabel('Nuclear Localization')
    %}
    
    %Plot NF and NMI for TM cells - MSN4
    
    figure(5)
    perm = [3,4,7,8] ;
    legend_vec = legend_vec_YFP(perm);
    cmap = cmap_YFP(perm,:);
    plot_params = {'linewidth',1.5,'LineStyle','-'};
    plt_grp = zeros(length(legend_vec),1);
    for jj = 1:length(perm)
        all_tracks = all_tracks_vec{perm(jj)};
        all_times = all_times_vec{perm(jj)};
        color_val = cmap(jj,:);
        plt_grp(jj) = hggroup;
        for ph = 1: length(phases)
            tracks = all_tracks.(phases{ph});
            timevals = all_times.(phases{ph});
            p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
            set(p,'Parent',plt_grp(jj))
        end
    end
    
    
    hleg = legend(plt_grp,legend_vec) %,'Location','NE');
    htitle = get(hleg,'Title');
    %set(htitle,'String','Condition')
    title('SC Cells, YFP Channel')
    xlabel('time')
    ylabel('Nuclear Localization')
    %}
    
    figure(6)
    perm = [3,4,7,8] ;
    legend_vec = legend_vec_YFP(perm);
    cmap = cmap_YFP(perm,:);
    plot_params = {'linewidth',1.5,'LineStyle','-'};
    plt_grp = zeros(length(legend_vec),1);
    for jj = 1:length(perm)
        all_tracks = all_tracks_vec{perm(jj)};
        all_times = all_times_vec{perm(jj)};
        color_val = cmap(jj,:);
        plt_grp(jj) = hggroup;
        for ph = 1: length(phases)
            tracks = all_tracks.(phases{ph});
            timevals = all_times.(phases{ph});
            p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nmi','plot_params',plot_params);
            set(p,'Parent',plt_grp(jj))
        end
    end
    
    
    hleg = legend(plt_grp,legend_vec) %,'Location','NE');
    htitle = get(hleg,'Title');
    %set(htitle,'String','Condition')
    title('SC Cells, YFP Channel')
    xlabel('time')
    ylabel('Median Intensity')
    %}
    
    
    
    
    
end

return

% 
% %combine legend vectors
% legend_vec = [legend_vec_RFP;legend_vec_YFP];
% hleg = legend(plt_grp,legend_vec) %,'Location','NE');
% htitle = get(hleg,'Title');
% %set(htitle,'String','Condition')
% title('SC.Hog1, KL.Msn2, SC.Msn2 response in S.Cerevisiae cells')
% xlabel('time')
% ylabel('Nuclear Localization')
% %}

%% Just SC.MSN2

figure(2)
clf
hold on
plt_grp = []
channel = 'RFP'
legend_vec_RFP = {'(SC.MSN2) GD -> GD + 0.25M Sorb t6', 
'(SC.MSN2) SDC -> SDC + 0.25M Sorb t6', 
'(SC.MSN2) GD -> GD + 0.25M Sorb t12', 
'(SC.MSN2) SDC -> SDC + 0.25M Sorb t12', 
'(SC.MSN2) GD -> GD + 0.25M Sorb t18', 
'(SC.MSN2) SDC -> SDC + 0.25M Sorb t18', 
'(SC.MSN2) 2uM NMPP1 -> 0.25M Sorb t6 + 2uM NMPP1', 
'(SC.MSN2) 2uM NMPP1 -> 0.25M Sorb t18 + 2uM NMPP1'}
fname_save = '20150204_processed_data_SC.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap = [0,0,1;  %Blue
 1,0,0;  %Red
  0,0.5,1;  %Light blue 
 1,0.35,0; %light red
 0,1,1; %Cyan
 1,0.7,0; %orange
 0.2,0.2,0.2; %gray
 0,0,0   %black
 ]; 

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,3,5,2,4,6] ;
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
N_RFP = length(perm);
for jj = 1:N_RFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end


%plot w/ legend
legend_vec = [legend_vec_RFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SC.MSN2 Gluc Drop followed by 0.25M Sorb')
xlabel('time')
ylabel('Nuclear Localization')
%axis([0,90,1.5,6.5])
%}


figure(3)
clf
hold on
plt_grp = []
channel = 'YFP'
legend_vec_YFP = {'(KL.MSN2) GD -> GD + 0.25M Sorb t6', 
'(KL.MSN2) SDC -> SDC + 0.25M Sorb t6', 
'(KL.MSN2) GD -> GD + 0.25M Sorb t12', 
'(KL.MSN2) SDC -> SDC + 0.25M Sorb t12', 
'(KL.MSN2) GD -> GD + 0.25M Sorb t18', 
'(KL.MSN2) SDC -> SDC + 0.25M Sorb t18', 
'(KL.MSN2) 2uM NMPP1 -> 0.25M Sorb t6 + 2uM NMPP1', 
'(KL.MSN2) 2uM NMPP1 -> 0.25M Sorb t18 + 2uM NMPP1'}
fname_save = '20150204_processed_data_SC.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap = [0,0,1;  %Blue
 1,0,0;  %Red
  0,0.5,1;  %Light blue 
 1,0.35,0; %light red
 0,1,1; %Cyan
 1,0.7,0; %orange
 0.2,0.2,0.2; %gray
 0,0,0   %black
 ];  


%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,3,5,2,4,6] ;
legend_vec_YFP = legend_vec_YFP(perm);
cmap = cmap(perm,:);
plot_params = {'linewidth',1.5,'LineStyle',':'};
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

legend_vec = [legend_vec_RFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 Gluc Drop followed by 0.25M Sorb')
xlabel('time')
ylabel('Nuclear Localization')
%}


return

%% Mean traces 0.5% plots

figure(3)
clf
hold on
plt_grp = []
channel = 'RFP'
legend_vec_RFP = { 
'(SC.MSN2) 2% -> 2% + 0.5M Glu', 
'(SC.MSN2) 2% -> 2% + 0.5M Sorb', 
'(SC.MSN2) 0.5% -> 0.5% + 0.5M Glu', 
'(SC.MSN2) 0.5% -> 0.5% + 0.5M Sorb', 
'(SC.MSN2) no glu -> .5M glu', 
'(SC.MSN2) no glu -> 0.5M Sorb', 
'(SC.MSN2) 2% -> 0.5% Glu', 
'(SC.MSN2) 2% -> 0% Glu'}
fname_save = '20150204_processed_data_SC.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap = [1,0,0;  %Red
 0,0,1;  %Blue
 1,0,0;  %Red
 0,0,1;  %Blue
 1,0,0;  %Red
 0,0,1;  %Blue
 0.5,0.5,0.5; %gray
 0,0,0   %black
 ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [3,4,7,8] ;
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
N_RFP = length(perm);
for jj = 1:N_RFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

channel = 'YFP'
legend_vec_YFP = {
'(KL.MSN2) 2% -> 2% + 0.5M Glu', 
'(KL.MSN2) 2% -> 2% + 0.5M Sorb', 
'(KL.MSN2) 0.5% -> 0.5% + 0.5M Glu', 
'(KL.MSN2) 0.5% -> 0.5% + 0.5M Sorb', 
'(KL.MSN2) no glu -> .5M glu', 
'(KL.MSN2) no glu -> 0.5M Sorb', 
'(KL.MSN2) 2% -> 0.5% Glu', 
'(KL.MSN2) 2% -> 0% Glu'}
fname_save = '20150204_processed_data_SC.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap = [1,0,0;  %Red
 0,0,1;  %Blue
 1,0,0;  %Red
 0,0,1;  %Blue
 1,0,0;  %Red
 0,0,1;  %Blue
 0.5,0.5,0.5; %gray
 0,0,0   %black
 ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [3,4,7,8] ;
legend_vec_YFP = legend_vec_YFP(perm);
cmap = cmap(perm,:);
plot_params = {'linewidth',1.5,'LineStyle',':'};
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj+N_RFP) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj+N_RFP))
    end
end


%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SD 0.5% Glucose')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,90,1.5,6.5])
%}


%% Mean traces no gluc plots

figure(4)
clf
hold on
plt_grp = []
channel = 'RFP'
legend_vec_RFP = { 
'(SC.MSN2) 2% -> 2% + 0.5M Glu', 
'(SC.MSN2) 2% -> 2% + 0.5M Sorb', 
'(SC.MSN2) 0.5% -> 0.5% + 0.5M Glu', 
'(SC.MSN2) 0.5% -> 0.5% + 0.5M Sorb', 
'(SC.MSN2) no glu -> .5M glu', 
'(SC.MSN2) no glu -> 0.5M Sorb', 
'(SC.MSN2) 2% -> 0.5% Glu', 
'(SC.MSN2) 2% -> 0% Glu'}
fname_save = '20150204_processed_data_SC.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap = [1,0,0;  %Red
 0,0,1;  %Blue
 1,0,0;  %Red
 0,0,1;  %Blue
 1,0,0;  %Red
 0,0,1;  %Blue
 0.5,0.5,0.5; %gray
 0,0,0   %black
 ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [5,6,7,8] ;
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
N_RFP = length(perm);
for jj = 1:N_RFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

channel = 'YFP'
legend_vec_YFP = {
'(KL.MSN2) 2% -> 2% + 0.5M Glu', 
'(KL.MSN2) 2% -> 2% + 0.5M Sorb', 
'(KL.MSN2) 0.5% -> 0.5% + 0.5M Glu', 
'(KL.MSN2) 0.5% -> 0.5% + 0.5M Sorb', 
'(KL.MSN2) no glu -> .5M glu', 
'(KL.MSN2) no glu -> 0.5M Sorb', 
'(KL.MSN2) 2% -> 0.5% Glu', 
'(KL.MSN2) 2% -> 0% Glu'}
fname_save = '20150204_processed_data_SC.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap = [1,0,0;  %Red
 0,0,1;  %Blue
 1,0,0;  %Red
 0,0,1;  %Blue
 1,0,0;  %Red
 0,0,1;  %Blue
 0.5,0.5,0.5; %gray
 0,0,0   %black
 ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [5,6,7,8] ;
legend_vec_YFP = legend_vec_YFP(perm);
cmap = cmap(perm,:);
plot_params = {'linewidth',1.5,'LineStyle',':'};
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj+N_RFP) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj+N_RFP))
    end
end


%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SD no Glucose')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,90,1.5,6.5])
%}

return


%%

figure(2)
clf
hold on
 
channel = []
%legend_vec_RFP = {'(KL.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorb';'(KL.MSN2) t1: SDC, t2: 0.5M Sorb'}
legend_vec_Hog1 = {'Hog1-YFP 2% Gluc + .5M Gluc',
    'Hog1-YFP 2% Gluc + .5M Sorb',
    'Hog1-YFP .11M Sorb',
    'Hog1-YFP SDC'}
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
fname_save = '20150123_processed_data_SC_Hog1.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_Hog1 = all_times_vec;
%all_tracks_vec_Hog1 = all_tracks_vec;
%posvec_Hog1 = posvec;
cmap = [1,0,0;  %Red
 0,1,0;  %Green
 0,0,1;  %Blue
 0,0,0   %Black
 ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);

perm = [1,2,3,4] ;
legend_vec_Hog1 = legend_vec_Hog1(perm);

%Get normalization factor (mean of control during Post phase)

%index of normalizing condition
norm_ind = 4
%Phase for normalization
norm_ph = 'Post'
norm_val = get_normval(all_tracks_vec, all_times_vec, norm_ind, norm_ph, channel)
plot_params = {'linewidth',1.5,'LineStyle','--'}

for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','norm_val',norm_val,'plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

%add Legend vector
legend_vec = legend_vec_Hog1;
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SC.Hog1 response in S.Cerevisiae cells')
xlabel('time')
ylabel('Nuclear Localization')
%}



figure(3)
clf 
hold on

channel = 'RFP'
%legend_vec_RFP = {'(KL.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorb';'(KL.MSN2) t1: SDC, t2: 0.5M Sorb'}
legend_vec_Msn2_RFP = {'SC MSN2-RFP 2% Gluc + .5M Gluc',
    'SC MSN2-RFP 2% Gluc + .5M Sorb',
    'SC MSN2-RFP .11M Sorb',
    'SC MSN2-RFP SDC'}
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
fname_save = '20150123_processed_data_SC_Msn2.mat';
load([base_dir,fname_save],'all_timevals_vec','all_tracks_vec','posvec')
%all_times_vec_Msn2 = all_times_vec;
%all_tracks_vec_Msn2 = all_tracks_vec;
%posvec_Msn2 = posvec;
cmap = [1,0,0;  %Red
 0,1,0;  %Green
 0,0,1;  %Blue
 0,0,0   %Black
 ];  

%index of normalizing condition
norm_ind = 4
%Phase for normalization
norm_ph = 'Post'
norm_val = get_normval(all_tracks_vec, all_times_vec, norm_ind, norm_ph, channel)
plot_params = {'linewidth',1.5,'LineStyle','-'}

for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','norm_val',norm_val,'plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end


%add Legend vector
legend_vec = legend_vec_Msn2_RFP;
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SC.MSN2 response in S.Cerevisiae cells')
xlabel('time')
ylabel('Nuclear Localization')
%}


figure(4)
clf 
hold on

channel = 'YFP'
%legend_vec_RFP = {'(KL.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorb';'(KL.MSN2) t1: SDC, t2: 0.5M Sorb'}
legend_vec_Msn2_YFP = {'KL MSN2-YFP 2% Gluc + .5M Gluc',
    'KL MSN2-YFP 2% Gluc + .5M Sorb',
    'KL MSN2-YFP .11M Sorb',
    'KL MSN2-YFP SDC'}
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
fname_save = '20150123_processed_data_SC_Msn2.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_Msn2 = all_times_vec;
%all_tracks_vec_Msn2 = all_tracks_vec;
%posvec_Msn2 = posvec;
cmap = [1,0,0;  %Red
 0,1,0;  %Green
 0,0,1;  %Blue
 0,0,0   %Black
 ];  

%index of normalizing condition
norm_ind = 4
%Phase for normalization
norm_ph = 'Post'
norm_val = get_normval(all_tracks_vec, all_times_vec, norm_ind, norm_ph, channel)
plot_params = {'linewidth',1.5,'LineStyle',':'}

for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','norm_val',norm_val,'plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end


%combine legend vectors
legend_vec = legend_vec_Msn2_YFP;
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.Msn2 response in S.Cerevisiae cells')
xlabel('time')
ylabel('Nuclear Localization')
%}



return

%{
figure(1)
clf
hold on

wells_avail = wellvecSP.(species)
channels_to_show = {'RFP','YFP'}
wells_to_show.RFP = {'A8','C8','E8','B8','D8', 'F8'}
wells_to_show.YFP = {'A8','C8','E8','B8','D8', 'F8'}
legend_vec.RFP = {'SC.MSN2 t0: Gluc DO t6: Gluc DO + 0.25M Sorbitol',
    'SC.MSN2 t0: Gluc DO t12: Gluc DO + 0.25M Sorbitol',
    'SC.MSN2 t0: Gluc DO t18: Gluc DO + 0.25M Sorbitol',
    'SC.MSN2 t0: SDC t6: SDC + 0.25M Sorbitol', 
    'SC.MSN2 t0: SDC t12: SDC + 0.25M Sorbitol',
    'SC.MSN2 t0: SDC t18: SDC + 0.25M Sorbitol'} ;
    %'SC.MSN2 t0: 2uM 1-NM-PP1 t6: 2uM 1-NM-PP1 + 0.25M Sorbitol',
    %'SC.MSN2 t0: 2uM 1-NM-PP1 t18: 2uM 1-NM-PP1 + 0.25M Sorbitol'}
legend_vec.YFP = {'KL.MSN2 t0: Gluc DO t6: Gluc DO + 0.25M Sorbitol',
    'KL.MSN2 t0: Gluc DO t12: Gluc DO + 0.25M Sorbitol',
    'KL.MSN2 t0: Gluc DO t18: Gluc DO + 0.25M Sorbitol',
    'KL.MSN2 t0: SDC t6: SDC + 0.25M Sorbitol', 
    'KL.MSN2 t0: SDC t12: SDC + 0.25M Sorbitol',
    'KL.MSN2 t0: SDC t18: SDC + 0.25M Sorbitol'} ;
    %'SC.MSN2 t0: 2uM 1-NM-PP1 t6: 2uM 1-NM-PP1 + 0.25M Sorbitol',
    %'SC.MSN2 t0: 2uM 1-NM-PP1 t18: 2uM 1-NM-PP1 + 0.25M Sorbitol'}
%     'KL.MSN2 t0: 2uM 1-NM-PP1 t18: 2uM 1-NM-PP1 + 0.25M Sorbitol'};
% wells_to_show.RFP = {'A8','B8','C8','D8', 'E8','F8','G8','H8'}
% wells_to_show.YFP = {'A8','B8','C8','D8', 'E8','F8','G8','H8'}
% legend_vec.RFP = {'SC.MSN2 t0: Gluc DO t6: Gluc DO + 0.25M Sorbitol',
%     'SC.MSN2 t0: SDC t6: SDC + 0.25M Sorbitol', 
%     'SC.MSN2 t0: Gluc DO t12: Gluc DO + 0.25M Sorbitol',
%     'SC.MSN2 t0: SDC t12: SDC + 0.25M Sorbitol', 
%     'SC.MSN2 t0: Gluc DO t18: Gluc DO + 0.25M Sorbitol',
%     'SC.MSN2 t0: SDC t18: SDC + 0.25M Sorbitol', 
%     'SC.MSN2 t0: 2uM 1-NM-PP1 t6: 2uM 1-NM-PP1 + 0.25M Sorbitol',
%     'SC.MSN2 t0: 2uM 1-NM-PP1 t18: 2uM 1-NM-PP1 + 0.25M Sorbitol'}
% legend_vec.YFP = {'KL.MSN2 t0: Gluc DO t6: Gluc DO + 0.25M Sorbitol',
%     'KL.MSN2 t0: SDC t6: SDC + 0.25M Sorbitol', 
%     'KL.MSN2 t0: Gluc DO t12: Gluc DO + 0.25M Sorbitol',
%     'KL.MSN2 t0: SDC t12: SDC + 0.25M Sorbitol', 
%     'KL.MSN2 t0: Gluc DO t18: Gluc DO + 0.25M Sorbitol',
%     'KL.MSN2 t0: SDC t18: SDC + 0.25M Sorbitol', 
%     'KL.MSN2 t0: 2uM 1-NM-PP1 t6: 2uM 1-NM-PP1 + 0.25M Sorbitol',
%     'KL.MSN2 t0: 2uM 1-NM-PP1 t18: 2uM 1-NM-PP1 + 0.25M Sorbitol'};
%cmap.RFP = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
% 0,0,0;
% 1,0,0];
%cmap.RFP = jet(length(legend_vec.RFP));
cmap.RFP = [        0         0    0.6667
         0         0    1.0000
         0    0.3333    1.0000
    1.0000         0         0
         1.0000    0.3333         0
         1.0000    0.6667         0];
cmap.YFP = cmap.RFP;
lstyle.RFP = '-';
lstyle.YFP = ':';

kk = 0; %set plot group number
for ch = 1:length(channels_to_show)
    channel = channels_to_show{ch}
    wts = wells_to_show.(channel);
    perm = zeros(length(wts),1);
    for pp = 1:length(wts)
        well = wts{pp};
        perm(pp) = find(strcmp(wells_avail,well));
    end
    cmap_ch = cmap.(channel);
        
    for jj = 1:length(perm)
        all_tracks = all_tracks_vec{perm(jj)};
        all_times = all_times_vec{perm(jj)};
        color_val = cmap_ch(jj,:);
        kk = kk+1
        plt_grp(kk) = hggroup;
        for ph = 1: length(phases)
            tracks = all_tracks.(phases{ph});
            times = all_times.(phases{ph});
            p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5,'LineStyle',lstyle.(channel));
            set(p,'Parent',plt_grp(kk))
        end
    end
end

legend_vec = [legend_vec.RFP;legend_vec.YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 and SC.MSN2 response in S.Cerevisiae cells')
xlabel('time')
ylabel('Nuclear Localization')
%}




%{
figure(1)
clf
hold on
channel = 'RFP'
%legend_vec_RFP = {'42 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol',
%    '42 t1: SDC, t2: 0.5M Sorbitol', 
%    '42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}

legend_vec_RFP = {'KL.MSN2 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol','KL.MSN2 t1: SDC, t2: 0.5M Sorbitol'}

%fname_save = '20140703_processed_data_KL_RFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;
% cmap = [0,0,1;
% 0,1,0.5;
% 0,.3333,.8333;
% 0,0.6667,0.6667;
% 1,0,0;
% 0,0,0]
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
% 0,0,0;
% 1,0,0];
cmap = jet(length(legend_vec_RFP));
perm = 1:length(legend_vec_RFP);
%perm = [2,1,6,5];
for jj = 1:length(legend_vec_RFP)
all_tracks = all_tracks_vec{perm(jj)};
all_times = all_times_vec{perm(jj)};
color_val = cmap(jj,:);
plt_grp(jj) = hggroup;
for ph = 1: length(phases)
tracks = all_tracks.(phases{ph});
times = all_times.(phases{ph});
p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5);
set(p,'Parent',plt_grp(jj))
end
end
% legend_vec_YFP = {'(SC.MSN2) 39y 1 0.5M sorb';
%     '(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
%     '(SC.MSN2) 39y 2 0.5M sorb';
%     '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'};
channel = 'YFP'
%legend_vec_YFP = {'42 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol',
%    '42 t1: SDC, t2: 0.5M Sorbitol', 
%    '42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}

legend_vec_YFP = {'SC.MSN2 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol','SC.MSN2 t1: SDC, t2: 0.5M Sorbitol'}

cmap = jet(length(legend_vec_RFP));
perm = 1:length(legend_vec_RFP);
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
%0,0,0;
%1,0,0
% ];
%perm = [2,1];
%fname_save = '20140703_processed_data_KL_YFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
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
p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5,'LineStyle',':');
set(p,'Parent',plt_grp(kk))
end
end
%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 and SC.MSN2 response in SC cells')
xlabel('time')
ylabel('Nuclear Localization')

%}

%S.Cer double pertubation experiment for BBC poster DEC 2014
%{
figure(2)
clf
hold on
channel = 'RFP'
%legend_vec_RFP = {'(KL.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorb';'(KL.MSN2) t1: SDC, t2: 0.5M Sorb'}
legend_vec_RFP = {'(SC.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol',
    '(SC.MSN2) t1: SDC, t2: 0.5M Sorbitol'}
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
%fname_save = '20140703_processed_data_KL_RFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;
cmap = [0,0,1;  %Blue
 0,0,0    %Black
 ];  

%cmap = [1,0,0;  %Red
% 0,0,0;  %Black
% 0,0,1;  %Blue
% 0.3,0.5,0.5; %Grey
% ];
% 0,.3333,.8333;
% 0,0.6667,0.6667;
% 1,0,0;
% 0,0,0]
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
% 0,0,0;
% 1,0,0];
%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2] ;
legend_vec_RFP = legend_vec_RFP(perm)
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5);
        set(p,'Parent',plt_grp(jj))
    end
end

channel = 'YFP'
%legend_vec_YFP =  {'(SC.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorb';'(SC.MSN2) t1: SDC, t2: 0.5M Sorb'}
legend_vec_YFP = {'(KL.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol',
    '(KL.MSN2) t1: SDC, t2: 0.5M Sorbitol'} 
%    '(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '(KL.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
%cmap = [1,0,0;  %Red
% 0,0,0;  %Black
% 0,0,1;  %Blue
% 0.3,0.5,0.5; %Grey
% ];
cmap = [0,0,1;  %Blue
 0,0,0    %Black
 ]; 
perm = [1,2];
%perm = [1,2,3,5] ;
legend_vec_YFP = legend_vec_YFP(perm);
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
%0,0,0;
%1,0,0
% ];
%perm = [2,1];
%fname_save = '20140703_processed_data_KL_YFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_YFP = all_times_vec;
%all_tracks_vec_YFP = all_tracks_vec;
%posvec_YFP = posvec;
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    kk = jj+length(legend_vec_RFP); %step up plot group
    plt_grp(kk) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5,'LineStyle',':');
        set(p,'Parent',plt_grp(kk))
    end
end
%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 and SC.MSN2 response in S.Cerevisiae cells')
xlabel('time')
ylabel('Nuclear Localization')
%}

%{
%K.Lac double pertubation experiment for BBC poster DEC 2014
figure(2)
%Ensure K.Lac is selected as species
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'(KL.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorb';'(KL.MSN2) t1: SDC, t2: 0.5M Sorb'}
%legend_vec_RFP = {'(SC.MSN2) 42 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol',
%    '(SC.MSN2) 42 t1: SDC, t2: 0.5M Sorbitol', 
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
%fname_save = '20140703_processed_data_KL_RFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;
cmap = [0,0,1;  %Blue
 0,0,0    %Black
 ];  

%cmap = [1,0,0;  %Red
% 0,0,0;  %Black
% 0,0,1;  %Blue
% 0.3,0.5,0.5; %Grey
% ];
% 0,.3333,.8333;
% 0,0.6667,0.6667;
% 1,0,0;
% 0,0,0]
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
% 0,0,0;
% 1,0,0];
%cmap = jet(length(legend_vec_RFP));
perm = 1:length(legend_vec_RFP);
%perm = [1,2,3,5] ;
legend_vec_RFP = legend_vec_RFP(perm)
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5);
        set(p,'Parent',plt_grp(jj))
    end
end

channel = 'YFP'
legend_vec_YFP =  {'(SC.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorb';'(SC.MSN2) t1: SDC, t2: 0.5M Sorb'}
%legend_vec_YFP = {'(KL.MSN2) 42 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol',
%    '(KL.MSN2) 42 t1: SDC, t2: 0.5M Sorbitol', 
%    '(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '(KL.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
%cmap = [1,0,0;  %Red
% 0,0,0;  %Black
% 0,0,1;  %Blue
% 0.3,0.5,0.5; %Grey
% ];
cmap = [0,0,1;  %Blue
 0,0,0    %Black
 ]; 
perm = [1,2];
%perm = [1,2,3,5] ;
legend_vec_YFP = legend_vec_YFP(perm);
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
%0,0,0;
%1,0,0
% ];
%perm = [2,1];
%fname_save = '20140703_processed_data_KL_YFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_YFP = all_times_vec;
%all_tracks_vec_YFP = all_tracks_vec;
%posvec_YFP = posvec;
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    kk = jj+length(legend_vec_RFP); %step up plot group
    plt_grp(kk) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5,'LineStyle',':');
        set(p,'Parent',plt_grp(kk))
    end
end
%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 and SC.MSN2 response in KL cells')
xlabel('time')
ylabel('Nuclear Localization')

%}

%{

figure(3)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1';
'(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol';
'(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol';
'(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol';
'(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol';
'(SC.MSN2) 42 t1: 0.5M Sorbitol t2: 0.5M Sorbitol';
'(SC.MSN2) 42 t1: 0.5M Sorbitol t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
%fname_save = '20140703_processed_data_KL_RFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;
cmap = [%0,0,0;  %Black
 0,0,1;  %Blue
 0,0.5,1; %Lighter Blue
 0,.5,0;  %Green
 0,1,0;  %Light Green
 ];
% 0,.3333,.8333;
% 0,0.6667,0.6667;
% 1,0,0;
% 0,0,0]
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
% 0,0,0;
% 1,0,0];
%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [2,3,4,5] ;
legend_vec_RFP = legend_vec_RFP(perm)
for jj = 1:length(perm)
all_tracks = all_tracks_vec{perm(jj)};
all_times = all_times_vec{perm(jj)};
color_val = cmap(jj,:);
plt_grp(jj) = hggroup;
for ph = 1: length(phases)
tracks = all_tracks.(phases{ph});
times = all_times.(phases{ph});
p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5);
set(p,'Parent',plt_grp(jj))
end
end
% legend_vec_YFP = {'(SC.MSN2) 39y 1 0.5M sorb';
%     '(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
%     '(SC.MSN2) 39y 2 0.5M sorb';
%     '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'};
channel = 'YFP'
legend_vec_YFP = {'(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1';
'(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol';
'(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol';
'(KL.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol';
'(KL.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol';
'(KL.MSN2) 42 t1: 0.5M Sorbitol t2: 0.5M Sorbitol';
'(KL.MSN2) 42 t1: 0.5M Sorbitol t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
cmap = [%0,0,0;  %Black
 0,0,1;  %Blue
 0,0.5,1; %Lighter Blue
 0,.5,0;  %
 0,1,0;  %
 ];
perm = [2,3,4,5] ;
legend_vec_YFP = legend_vec_YFP(perm);
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
%0,0,0;
%1,0,0
% ];
%perm = [2,1];
%fname_save = '20140703_processed_data_KL_YFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_YFP = all_times_vec;
%all_tracks_vec_YFP = all_tracks_vec;
%posvec_YFP = posvec;
for jj = 1:length(perm)
all_tracks = all_tracks_vec{perm(jj)};
all_times = all_times_vec{perm(jj)};
color_val = cmap(jj,:);
kk = jj+length(legend_vec_RFP); %step up plot group
plt_grp(kk) = hggroup;
for ph = 1: length(phases)
tracks = all_tracks.(phases{ph});
times = all_times.(phases{ph});
p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5,'LineStyle',':');
set(p,'Parent',plt_grp(kk))
end
end
%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 and SC.MSN2 response in SC cells')
xlabel('time')
ylabel('Nuclear Localization')


figure(4)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1';
'(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol';
'(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol';
'(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol';
'(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol';
'(SC.MSN2) 42 t1: 0.5M Sorbitol t2: 0.5M Sorbitol';
'(SC.MSN2) 42 t1: 0.5M Sorbitol t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
%fname_save = '20140703_processed_data_KL_RFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;
cmap = [%0,0,0;  %Black
 0,0,1;  %Blue
 1,0,0; %Red
 1,.5,.25; %Orange
 %0,1,0;  %Light Green
 ];
perm = [2,6,7] ;
% 0,.3333,.8333;
% 0,0.6667,0.6667;
% 1,0,0;
% 0,0,0]
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
% 0,0,0;
% 1,0,0];
%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
legend_vec_RFP = legend_vec_RFP(perm)
for jj = 1:length(perm)
all_tracks = all_tracks_vec{perm(jj)};
all_times = all_times_vec{perm(jj)};
color_val = cmap(jj,:);
plt_grp(jj) = hggroup;
for ph = 1: length(phases)
tracks = all_tracks.(phases{ph});
times = all_times.(phases{ph});
p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5);
set(p,'Parent',plt_grp(jj))
end
end
% legend_vec_YFP = {'(SC.MSN2) 39y 1 0.5M sorb';
%     '(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
%     '(SC.MSN2) 39y 2 0.5M sorb';
%     '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'};
channel = 'YFP'
legend_vec_YFP = {'(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1';
'(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol';
'(KL.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol';
'(KL.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol';
'(KL.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol';
'(KL.MSN2) 42 t1: 0.5M Sorbitol t2: 0.5M Sorbitol';
'(KL.MSN2) 42 t1: 0.5M Sorbitol t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
cmap = [%0,0,0;  %Black
 0,0,1;  %Blue
 1,0,0; %Red
 1,.5,.25; %Orange
 %0,1,0;  %Light Green
 ];
perm = [2,6,7] ;
legend_vec_YFP = legend_vec_YFP(perm);
%cmap = [0,0,1;
% .5,0.25,0;
%0,.3333,.8333;
%0,0.6667,0.6667;
%0,0,0;
%1,0,0
% ];
%perm = [2,1];
%fname_save = '20140703_processed_data_KL_YFP.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
%all_times_vec_YFP = all_times_vec;
%all_tracks_vec_YFP = all_tracks_vec;
%posvec_YFP = posvec;
for jj = 1:length(perm)
all_tracks = all_tracks_vec{perm(jj)};
all_times = all_times_vec{perm(jj)};
color_val = cmap(jj,:);
kk = jj+length(legend_vec_RFP); %step up plot group
plt_grp(kk) = hggroup;
for ph = 1: length(phases)
tracks = all_tracks.(phases{ph});
times = all_times.(phases{ph});
p = plot_meanvalues(times,tracks,channel,color_val,0,'nf','linewidth',1.5,'LineStyle',':');
set(p,'Parent',plt_grp(kk))
end
end
%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 and SC.MSN2 response in SC cells')
xlabel('time')
ylabel('Nuclear Localization')






%{

%the site I want to image is B11 which is the first site in the stored
%data. A11 has a bad first image.

sites ={'A1','B1','C1','D1', 'E1', 'F1', 'G1', 'H1'};
wellvecSP.SC = {'A1','B1','C1','D1', 'E1', 'F1', 'G1', 'H1'};

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

%}


%{
%Normalize by norm_val
%plot
phase = 'Post'
site = 'A8'
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
title(['SC: SC.MSN2 nuclear localization. ', condition, '.'])
ylabel('Cells')
xlabel('Time')
set(gca,'XTick',1:length([sing_cell_tracks.(site).(phase).times])+0.5)
set(gca,'XTickLabel',sprintf('%0.0f|',[sing_cell_tracks.(site).(phase).times]))
caxis([-0.4,4])

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
title(['SC: KL.MSN2 nuclear localization. ', condition, '.'])
ylabel('Cells')
xlabel('Time')
set(gca,'XTick',1:length([sing_cell_tracks.(site).(phase).times])+0.5)
set(gca,'XTickLabel',sprintf('%0.0f|',[sing_cell_tracks.(site).(phase).times]))
caxis([-0.4,4])

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
xlabel('SC.MSN2(RFP)')
ylabel('KL.MSN2(YFP)')
title(['SC: ',condition])
axis([1,14,1,9])
%scatterplot of data altering color

figure(4)
plot([sing_cell_tracks.(site).(phase).times],corr_vec,'LineWidth',3)
title(['SC: ',condition,'. Correlation SC.MSN2 n.loc. to KL.MSN2 n.loc'])
xlabel('Time')
ylabel('Corr')
axis([10,70,0,1])

%set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
%cmap = jet(length(legend_vec));
cmap = [1,0,0;
0,0,0];

%}


%{
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

%}

%}


