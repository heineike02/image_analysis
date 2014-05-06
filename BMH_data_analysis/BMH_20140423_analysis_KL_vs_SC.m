function all_tracks = BMH_20140423_analysis_KL_vs_SC()
%Separate plotting from data collection
%Save data in between data collection
%Plot median intensity
%Plot number of cells
%Spot cells with bright field
%Collect Cell Size



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

%line: 37 KL.BGLII , doubledash, 38: KL.COXI, dots: 39: ag.TEF1

%input information from experiment here
species_cell = {'SC'} %{'SC'}  %
phases =  {'Pre','Post'} %,'Post'} 
shift_spacing = [0,4]

channels = {'BF','RFP','YFP'}
channel_to_image = 'RFP'

dt = 4
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

bgimg = 0; 
%1 if you have a background image, 0 if you don't. 
%Gets background image depending on channel to image
%Collect background images using micromanager 
if bgimg == 0
    imbg = 1 % default
else
    %imbg = imread([base_dir,'BG\',channel_to_image,'_p1.tiff']);
    if strcmp(channel_to_image,'RFP')
         imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140423\RFP_BG\img_000000000_Default_000.tiff');
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

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc
posvec = {'p1','p2','p3';
'p4','p5','p6';
'p7','p8','p9';
'p10','p11','p12';
'p13','p14','p15';
'p16','p17','p18';
'p19','p20','p21';
'p22','p23','p24'};

%set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
cmap = cool(length(dosevec));

figure(1)
clf 
hold on
for jj = 1:length(dosevec)
    pos_fnamesSP.SC = posvec(jj,:)   %{'p16','p17','p18'} %, 'p14','p15'}
    [all_tracks,all_times] = KL_vs_SC_analysis(ipdir,storeim,fname_conv,op_amp,std_threshSP,species_cell,phases, shift_spacing,imdirPhase, maxdisp_1x,pos_fnamesSP,channels,channel_to_image,dt,imbg)
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0);
        set(p,'Parent',plt_grp(jj))
    end
end

%convert dosevec to string
for jj = 1: length(dosevec)
    dosevec_str{jj} = num2str(dosevec(jj),'%0.2f')
end
hleg = legend(plt_grp,dosevec_str,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Tm (ug/ml)')

title('MSN2 Nuclear Localization after Tm treatment')




end



