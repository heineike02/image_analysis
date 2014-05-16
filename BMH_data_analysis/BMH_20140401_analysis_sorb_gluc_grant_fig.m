function BMH_20140401_analysis_sorb_gluc_grant_fig()
% Glucose Dropout and Osmotic stress in Lactis and Cerevisiae
% Visualize individual and mean nuclear localization each condition
% in S.C. and K. Lactis

% Used 1.5x optical zoom.  
% Need to have argument for channel to measure / output
% Change jacobs image code output to sort in image order. 
% Too many tracks - something wierd going on

% Had to use ".tiff" for several images
% Note: to make training image use ginput function

%Plot number of cells
%Spot cells with bright field
%Collect Cell Size


ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(path,'..\..\image_analysis')

% Filename convention
% 'JSO' or 'Micromanager'
% For example of Micromanager filename convention see BMH_20140127 analysis
% files
fname_conv = 'JSO'

%imdirPhase.Pre = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_pre_1\'
%imdirPhase.Post = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_post_1\'

base_dir = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140401\'
imdirPhase.Pre = [base_dir,'OS_GD_Pre\']
imdirPhase.Post = [base_dir,'OS_GD_Post\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%line: 37 KL.BGLII , doubledash, 38: KL.COXI, dots: 39: ag.TEF1

%input information from experiment here
species_cell = {'SC'} %{'SC'}  %Right now not properly cycling through each species
phases =  {'Pre','Post'} %,'Post'} 
shift_spacing = [0,4]    

channels = {'BF','RFP','BF'}
channel_to_image = 'RFP'

dt = 2.5
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
if bgimg == 0
    imbg = 1 % default
else
    %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
    if strcmp(channel_to_image,'RFP')
        imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140401\OS_GD_post\RFP_p5_t2.tiff');
    elseif strcmp(channel_to_image,'YFP')
        imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140401\OS_GD_post\YFP_BG.tif');
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


%Pos 1-3: SC GD
%Pos 4-6: KL GD
%Pos 7-9: SC Sorb
%Pos 10-12: KL Sorb
%Pos 13-15: SC  NaCL
%Pos 16-18: KL NaCL
%Pos 19-21: SC Cont
%Pos 22-24: KL Cont

cond_vec = {'Glucose Dropout 0.2%','Sorbitol 0.25M','Control'};

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc
posvec.SC = {'p1','p2','p3';
'p7','p8','p9';
'p19','p20','p21'};

posvec.KL = {'p4','p5','p6';
'p10','p11','p12';
'p22','p23','p24'};

%Obtain and store data for each dose (note: only need to do this once) 

get_data = 0

all_tracks_vec = [];
all_times_vec = [];
if get_data == 1 
    for jj = 1:length(cond_vec)
       pos_fnamesSP.SC = posvec.SC(jj,:)   %{'p16','p17','p18'} %, 'p14','p15'}
       pos_fnamesSP.KL = posvec.KL(jj,:)
       [all_tracks,all_times] = KL_vs_SC_analysis(ipdir,storeim,fname_conv,op_amp,std_threshSP,species_cell,phases, shift_spacing,imdirPhase, maxdisp_1x,pos_fnamesSP,channels,channel_to_image,dt,imbg)
       all_tracks_vec{jj} = all_tracks;
       all_times_vec{jj} = all_times;
    end
    %store data
    save([base_dir,'grant_prop_processed_data_SC.mat'],'all_times_vec','all_tracks_vec','cond_vec','posvec')
else
    %retreive data
    load([base_dir,'grant_prop_processed_data_SC.mat'],'all_times_vec','all_tracks_vec','cond_vec','posvec')   
end

%set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
cmap = [[0,0,1];[1,0,0];[0,0,0]];

figure(1)
clf 
hold on
for jj = 1:length(cond_vec)
    all_tracks = all_tracks_vec{jj};
    all_times = all_times_vec{jj};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,1,'nf');
        set(p,'Parent',plt_grp(jj))
    end
end



hleg = legend(plt_grp,cond_vec) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Condition')
title('S.Cerevisae MSN2 Nuclear Localization','FontSize',12,'FontWeight','demi')
%axis([0,65,1.25,3.5]) %K. Lactis
axis([0,65,1.25,6.5]) %S. Cerivisiae

figure(2)
clf 
hold on
for jj = 1:length(cond_vec)
    all_tracks = all_tracks_vec{jj};
    all_times = all_times_vec{jj};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,1,'nmi');
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,cond_vec) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Sorbitol (M)')
title('MSN2 Nuclear Localization after sorbitol treatment')


end


