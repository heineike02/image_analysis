function all_tracks = BMH_20140516_analysis_gluc_dose_resp()

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
% 'JSO','Micromanager' or 'HTC_Nikon'
% For example of Micromanager filename convention see BMH_20140127 analysis
% files
fname_conv = 'HTC_Nikon'

%imdirPhase.Pre = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_pre_1\'
%imdirPhase.Post = 'C:\Users\Ben\Documents\Data\PKA_Project\20140123\10_27_28_GD_GA_post_1\'

base_dir = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140516_Gluc_drop_dose_resp\'
imdirPhase.Pre = [base_dir,'gluc_drop_pre\']
imdirPhase.Post = [base_dir,'gluc_drop_post\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%input information from experiment here
species = 'SC' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species

channels = {'BF','RFP_cube'}
channel_to_image = 'RFP_cube'

phases =  {'Pre','Post'} %,'Post'} 
shift_spacing = [0,12]    
%timestep method
%time_calc:  Tells the program how to calculate each time value.  If the
%input is a number then that is interpreted as a dt between images
%If the input is a filename, then that is interpreted as metadata that gives
%exact time values.  

%time_calc_phase.Pre = 4
%time_calc_phase.Post = 4
time_calc_phase.Pre = [imdirPhase.Pre, 'acqdat.txt']
time_calc_phase.Post = [imdirPhase.Post, 'acqdat.txt']


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
    if strcmp(channel_to_image,'RFP_cube')
         imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140430\BG\RFP_cube_BG\img_000000000_Default_000.tif');
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


%Pos 1-3: SC 0.5 M Sorb
%Pos 4-5: KL 0.5 M Sorb
%Pos 6-8: SC 0.25 M Sorb
%Pos 9-10: KL 0.25 M Sorb
%Pos 11-13: SC 0.125 M Sorb
%Pos 14-15: KL 0.125 M Sorb
%Pos 16-18: SC 0.0625 M Sorb
%Pos 19-20: KL 0.0625 M Sorb
%Pos 21-23: SC 0.03125 M Sorb
%Pos 24-25: KL 0.03125 M Sorb
%Pos 26-28: SC Cont
%Pos 29-30: KL Cont

dosevec = [0.5,0.25,0.125,0.0625,0.03125,0];

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc
posvec.SC = {'p1','p2','p3';
'p6','p7','p8';
'p11','p12','p13';
'p16','p17','p18';
'p21','p22','p23';
'p26','p27','p28'};

posvec.KL = {'p4','p5';
'p9','p10';
'p14','p15';
'p19','p20';
'p24','p25';
'p29','p30'};

%Obtain and store data for each dose (note: only need to do this once) 

get_data = 1;

all_tracks_vec = [];
all_times_vec = [];
if get_data == 1 
    for jj = 1:length(dosevec)
       for ph = 1:length(phases);
            phase = phases{ph};
            imdir = imdirPhase.(phase);
            shift_spacing_ph = shift_spacing(ph);
            time_calc = time_calc_phase.(phase);
            pos_fnamesSP.SC = posvec.SC(jj,:);   %{'p16','p17','p18'} %, 'p14','p15'}
            pos_fnamesSP.KL = posvec.KL(jj,:);
            [tracks,times] = KL_vs_SC_analysis(ipdir,storeim,fname_conv,op_amp,std_threshSP,species, imdir, maxdisp_1x,pos_fnamesSP,channels,channel_to_image,time_calc,imbg);
            all_tracks.(phase) = tracks;
            all_times.(phase) = times + shift_spacing(ph);
       end
       all_tracks_vec{jj} = all_tracks;
       all_times_vec{jj} = all_times;
    end
    %store data
    save([base_dir,'sorb_doseresp_processed_data_SC.mat'],'all_times_vec','all_tracks_vec','dosevec','posvec')
else
    %retreive data
    load([base_dir,'sorb_doseresp_processed_data_SC.mat'],'all_times_vec','all_tracks_vec','dosevec','posvec')   
end

%set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
cmap = cool(length(dosevec));

%convert dosevec to string
for jj = 1: length(dosevec)
    dosevec_str{jj} = num2str(dosevec(jj),'%0.2f');
end

figure(1)
clf 
hold on
for jj = 1:length(dosevec)
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



hleg = legend(plt_grp,dosevec_str) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Sorbitol (M)')
title('SC: MSN2 Nuclear Localization after sorbitol treatment')


figure(2)
clf 
hold on
for jj = 1:length(dosevec)
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

hleg = legend(plt_grp,dosevec_str) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Sorbitol (M)')
title('SC MSN2 median intensity')


figure(3)
clf 
hold on
for jj = 1:length(dosevec)
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

hleg = legend(plt_grp,dosevec_str); %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Sorbitol (M)')
title('SC: Number of cells identified')


end



