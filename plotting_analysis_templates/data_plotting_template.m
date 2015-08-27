function all_tracks = data_plotting_template()
%

profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150608_deltaT_GD_first\'
%input information from experiment here
species = 'SC' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species

phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0,12.5]    

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


%Strain BMH42 - SC.MSN2-RFP,  KL.MSN2-YFP 
%wellvecSP.SC = {'A9','B9','C9','D9','E9','F9','G9','H9'};
%

%all locations had 4 sites
% Glucose Dropout (0.1 % from 2%) followed by Osmo shock (0.35M) at various times.
% osmo balanced with sorbitol relative to osmotic stress condition.
%
%
%A8	t9: GD + sorb   
%B8	t9: GD      
%C8	t9: SDC + sorb
%D8	t9: GD  	t15: sorb
%E8	t9: GD      t21: sorb
%F8	t9: GD      t33: sorb
%G8	t9: GD      t45: sorb
%H8	t9: GD      t69: sorb



%% SC.Msn2-RFP Mean Time Traces
figure(1)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'t9: GD + sorb',
't9: GD'}

fname_save = 'example_processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [ 0,0,1;  %Blue
    1,0,0; %Red
    ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2]  ;
N_RFP = length(perm)
legend_vec_RFP_plot = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_RFP_plot),1);
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

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('GD and Osmo shock')
xlabel('time')
ylabel('Nuclear Localization')

%% KL.MSN2-YFP  Mean Time Traces
figure(2)
clf 
hold on
channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP
%perm same as above
legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_YFP_plot),1);
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


hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('GD and Osmo shock, varying times')
xlabel('time')
ylabel('Nuclear Localization')



%% SC.MSN2-RFP plot individual cells 
figure(3)
clf
hold on
channel = 'RFP'
title_vec = {'t9: GD + sorb',
't9: GD'};

%already loaded

perm = [1,2] ;
N_RFP = length(perm);
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);

for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    subplot(1,2,jj)
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_individual_values(timevals,tracks,channel,'nf');
    end
    axis([0,100,0,10])
    title(title_vec{jj})
end


