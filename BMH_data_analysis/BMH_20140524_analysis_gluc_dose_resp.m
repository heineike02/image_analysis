function all_tracks = BMH_20140524_analysis_gluc_dose_resp()

%not all wells have the same number of time points - the scope failed 

%to import an image sequence in imagej use regular expressions
%eg. WellA05_Seq0000t[0-9]*xy1c1.*  WellD05_Seq[0-9]*t[0-9]*xy1c1.*
%c1 is RFP


%Spot cells with bright field
%Collect Cell Size

% Change jacobs image code output to sort in image order. 

ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(path,'..\..\image_analysis')

% Filename convention
% 'JSO','Micromanager', or 'HCS_Nikon'
% For example of Micromanager filename convention see BMH_20140127 analysis
% files
fname_conv = 'Micromanager'

%imdirPhase.Pre = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_pre_1\'
%imdirPhase.Post = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_post_1\'

base_dir = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140524_gluc_drop\'
imdirPhase.Pre = [base_dir,'Pre\']
imdirPhase.Post = [base_dir,'Post\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%input information from experiment here
species = 'SC' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species
species_cell = {'SC','KL'}

fname_saveSP.KL = 'gluc_drop_doseresp_procdata_KL.mat';
fname_saveSP.SC = 'gluc_drop_doseresp_procdata_SC.mat';

channels = {'BF','RFP'}
channel_to_image = 'RFP'

phases =  {'Pre','Post'} %,'Post'} 
shift_spacing = [0,14]    
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
addpath('C:/Users/Ben/Documents/GitHub/image_analysis/jsonlab');
time_calc_phase.Pre  = 'metadata.txt'
time_calc_phase.Post = 'metadata.txt'

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
    imbg = 1 % default
else
    %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
    if strcmp(channel_to_image,'RFP')
         imbg = imread([base_dir,'\BG_RFP_gluc_drop\img_000000000_RFP_000.tif']);
         imbg = imbg'; %Micromanager images are the transpose of images taken by the image collection scripts.
    elseif strcmp(channel_to_image,'YFP')
         %imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\YFP_BG\img_000000000_Default_000.tiff');
    else
         'Error: imbg not assigned'
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


%A9 sites 1-4 SC 2% to 0.111M Sorb
%A10 sites 1-4 KL 2% to 0.111M Sorb
%B9 sites 1-4 SC 2% to no gluc
%B10 sites 1-4 KL 2% to no gluc
%C9 sites 1-4 SC 2% to 0.006M + 2% gly
%C10 sites 1-4 KL 2% to 0.006M + 2% gly
%D9 sites 1-4 SC 2% to 0.006M + 0.111M gly
%D10 sites 1-4 KL 2% to 0.006M + 0.111M gly
%E9 sites 1-4 SC 2% to 0.006M + Sorb
%E10 sites 1-4 KL 2% to 0.006M + Sorb
%F9 sites 1-4 SC 2% to 0.006M gluc
%F10 sites 1-4 KL 2% to  0.006M gluc 
%G9 sites 1-4 SC 2% to 0.25M Sorb + 2%gluc
%G10 sites 1-4 KL 2% to 0.25M Sorb + 2%gluc
%H9 sites 1-4 SC 2% to 0.111M
%H10 sites 1-4 KL 2% to 0.111M

legend_vec = {'0.111M Sorb','No Gluc' ,'0.006M Gluc + 2% Gly', '0.006M Gluc + 0.105M Gly', '0.006M Gluc + 0.105M Sorb', '0.006M Gluc', '0.25M Sorb + 2% Gluc', 'SDC 2%'}

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc

wellvecSP.SC = {'A9','B9','C9','D9','E9','F9','G9','H9'};
wellvecSP.KL = {'A10','B10','C10','D10','E10','F10','G10','H10'}
Nsites = 4;

for sp = 1:length(species_cell);
    wellvec = wellvecSP.(species_cell{sp});
    for jj = 1:length(wellvec);
        for kk = 1:Nsites;
            posvecSP.(species_cell{sp}){jj,kk} = [wellvec{jj},'-Site_',num2str(kk-1)];
        end
    end
end

%For some reason A7-Site_2 is bad, remove it from the list
%posvecSP.SC{1,3} = 'NA'

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

get_data = 1;

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
            all_times.(phase) = times + shift_spacing(ph);
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
        p = plot_meanvalues(times,tracks,color_val,0,'nf');
        set(p,'Parent',plt_grp(jj))
    end
end



hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Sorbitol (M)')
title('SC MSN2 Nuclear Localization after sorbitol treatment')


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
        p = plot_meanvalues(times,tracks,color_val,0,'nmi');
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Sorbitol (M)')
title('SC MSN2 median intensity')


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
set(htitle,'String','Sorbitol (M)')
title('SC Number of cells identified')


end



