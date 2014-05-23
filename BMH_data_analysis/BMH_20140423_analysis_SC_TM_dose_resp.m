function all_tracks = BMH_20140423_analysis_SC_TM_dose_resp()
%pos25 - difficult combinatorics - check it.  
%Also didn't use background image - should run again to improve data.
%Jacobs images are a reflection of the way micromanager saves images.

% Need to have argument for channel to measure / output
% Change jacobs image code output to sort in image order. 
% Too many tracks - something wierd going on

% Had to use ".tiff" for several images
% Note: to make training image use ginput function

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

base_dir = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\'
imdirPhase.Pre = [base_dir,'TM_doseresp_pre\']
imdirPhase.Post = [base_dir,'TM_doseresp_post\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']


%input information from experiment here
species = 'SC' %{'SC'}  %

channels = {'BF','RFP','YFP'}
channel_to_image = 'RFP'

phases =  {'Pre'} %,'Post'} %,'Post'}  
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
op_amp =  '1x'
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
         imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\RFP_BG\img_000000000_Default_000.tif');
         imbg = imbg';
    elseif strcmp(channel_to_image,'YFP')
         imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\YFP_BG\img_000000000_Default_000.tiff');
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


%Pos 1-3: 0 (control)
%Pos 4-6: 0.08
%Pos 7-9: 0.16
%Pos 10-12: 0.3
%Pos 13-15: 0.625
%Pos 16-18: 1.25
%Pos 19-21: 2.5
%Pos 22-24: 5ug/ml TM
%Pos 25-27: HSP12-YFP
%Pos 28-30: UPRE-YFP
%Pos 31-33: 11-38 NACL 0.25M

dosevec = [0.00,0.078125,0.15625,0.3125,0.625,1.25,2.50,5.00];
%convert dosevec to string
for jj = 1: length(dosevec)
    dosevec_str{jj} = num2str(dosevec(jj),'%0.2f');
end
legend_str = [dosevec_str,'HSP12','UPRE4x','NaCl 0.25M'];

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc
posvec = {'p1','p2','p3';
'p4','p5','p6';
'p7','p8','p9';
'p10','p11','p12';
'p13','p14','p15';
'p16','p17','p18';
'p19','p20','p21';
'p22','p23','p24';
'p25','p26','p27';
'p28','p29','p30';
'p31','p32','p33'};

%Obtain and store data for each dose (note: only need to do this once) 

get_data = 0

all_tracks_vec = [];
all_times_vec = [];
if get_data == 1 
    for jj = 1:length(posvec)
       if (jj==9)||(jj==10)  %UPRE and HSP12 controls
          channel_to_image = 'YFP' 
       else 
          channel_to_image = 'RFP'
       end
       
       for ph = 1:length(phases);
           phase = phases{ph};
           imdir = imdirPhase.(phase);
           shift_spacing_ph = shift_spacing(ph);
           time_calc = time_calc_phase.(phase);
           pos_fnamesSP.SC = posvec(jj,:)   %{'p16','p17','p18'} %, 'p14','p15'}
           [tracks,times] = KL_vs_SC_analysis(ipdir,storeim,fname_conv,op_amp,std_threshSP,species,imdir, maxdisp_1x,pos_fnamesSP,channels,channel_to_image,time_calc,imbg);
           all_tracks.(phase) = tracks;
           all_times.(phase) = times + shift_spacing(ph);
       end
       all_tracks_vec{jj} = all_tracks;
       all_times_vec{jj} = all_times;
    end
    %store data
    save([base_dir,'Pre_test.mat'],'all_times_vec','all_tracks_vec','dosevec','posvec')
else
    %retreive data
    load([base_dir,'Pre_test.mat'],'all_times_vec','all_tracks_vec','dosevec','posvec')   
end

%set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
cmap = cool(length(dosevec));
cmap = [cmap;[1,0,0];[1,0.5,0];[0,0,0]]; %add colors for controls

figure(1)
clf 
hold on
for jj = 1:length(posvec)
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



hleg = legend(plt_grp,legend_str); %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Tm (ug/ml)')
title('MSN2 Nuclear Localization after Tm treatment')


figure(2)
clf 
hold on
for jj = 1:length(posvec)
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

hleg = legend(plt_grp,legend_str); %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Tm (ug/ml)')
title('Median intensity after Tm treatment')

figure(3)
clf 
hold on
for jj = 1:length(posvec)
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

hleg = legend(plt_grp,legend_str); %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Tm (ug/ml)')
title('Number of cells identified')

end



