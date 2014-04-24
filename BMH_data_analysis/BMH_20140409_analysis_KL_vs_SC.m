function all_tracks = BMH_20140409_analysis_KL_vs_SC()

% Used 1.5x optical zoom.  
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

base_dir = 'C:\Users\Ben\Box Sync\Data\PKA_Project\20140409\'
imdirPhase.Pre = [base_dir,'TM_Sorb_H2O2_Pre\']
imdirPhase.Post = [base_dir,'TM_Sorb_H2O2_Post\']
%imdirPhase.Rep = [base_dir,'Gluc_Rep_Post\']

%input information from experiment here
species_cell = {'SC'} %{'SC'}  %
phases =  {'Pre','Post'} 
shift_spacing = [0,4]

%Pos 1-2: SC Cont 
%Pos 3-4: KL YFP Cont 
%Pos 5-6: KL RFP Cont
%Pos 7-8: SC TM 5ug/ml
%Pos 9-10: KL YFP TM 5ug/ml
%Pos 11-12: KL RFP TM 5ug/ml 
%Pos 13-14: SC SOrb 0.25M
%Pos 15-16: KL YFP Sorb
%Pos 17-18: KL RFP Sorb
%Pos 19-20: SC H2O2 1.0mM
%Pos 21-22: KL YFP H2O2
%Pos 23-24: KL RFP H2O2

%Micromanager positions: Pos0, Pos1, etc.  
%JSO positions p1,p2,etc
%pos_fnamesSP.SC = {'Pos0','Pos1','Pos2','Pos3','Pos4','Pos5','Pos6'};
%pos_fnamesSP.KL = {'Pos3','Pos4','Pos5'};
pos_fnamesSP.SC = {'p17','p16'} %, 'p14','p15'}
pos_fnamesSP.KL = {'p23','p24'} 
channels = {'YFP','RFP','BF'}
channel_to_image = 'RFP'
dt = 2.5

summary_title = 'KL H2O2 1.0(?)mM (blue), Sorb 0.25M (red), Control(black)'
acolor_summary = 'b'

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
         imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140409\TM_Sorb_H2O2_post\RFP_BG.tif');
    elseif strcmp(channel_to_image,'YFP')
         imbg = imread('C:\Users\Ben\Box Sync\Data\PKA_Project\20140409\TM_Sorb_H2O2_post\YFP_BG.tif');
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

all_tracks = KL_vs_SC_analysis(ipdir,storeim,fname_conv,op_amp,std_threshSP,species_cell,phases, shift_spacing,imdirPhase, maxdisp_1x,pos_fnamesSP,channels,channel_to_image,dt,imbg,summary_title,acolor_summary)

end



