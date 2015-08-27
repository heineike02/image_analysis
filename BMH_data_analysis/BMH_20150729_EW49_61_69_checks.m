function [] = BMH_20150729_EW49_61_69_checks()
%

profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150729_EW049_61_69_check_osmo\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    


%Strain BMH42 - SC.MSN2-RFP,  KL.MSN2-YFP 
%

%Exp 1
%all locations had 4 sites
% Osmo shock 0.5M Sorbitol added after Pre phase (approx 10 min)
% Row #     Strn GD         Osmo        Fluor
% GD / Osmo
% A8 / A8	49 1 GD	        0.5M final	rfp
% B8 / B8	49 2 GD	        0.5M final	rfp
% C8 / C8	37   GD	        0.5M final	both
% D8 / D8	37   GD 0.1%	0.25M final	both
% E8 / E8	61-1 GD	        0.5M final	yfp
% F8 / F8	61-1 GD 0.1%	0.25M final	yfp
% G8 /G8	69-1 GD	        0.5M final	yfp
% H8 / H8 	69-1 GD 0.1%	0.25M final	yfp




%% Osmo shock plots

%RFP channel 

figure(1)
clf
hold on

channel = 'RFP'
legend_vec_RFP = {'BMH37 0.5M','BMH37 0.25M','EW49 1', 'EW49 2'}
cmap_RFP = [ 0,0,0;  %Black
    0.5,0.5,0.5; %Grey
    1,0,0;  %Red
    0.5,0,0 %Red 2
        ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3,4]  ;
N_RFP = length(perm)  
    
fname_save = 'processed_data_Osmo_37.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.BMH37 = all_tracks_vec;
all_times_vec_comb.BMH37 = all_times_vec;

fname_save = 'processed_data_Osmo_49.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW49 = all_tracks_vec;
all_times_vec_comb.EW49 = all_times_vec;

data_file_list = {'BMH37','EW49'};
clear('all_tracks_vec');
clear('all_times_vec');
all_tracks_vec = {};
all_times_vec = {};
for jj = 1:length(data_file_list);
    all_tracks_vec = [all_tracks_vec , all_tracks_vec_comb.(data_file_list{jj})];
    all_times_vec = [all_times_vec , all_times_vec_comb.(data_file_list{jj})];
end



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
title('Osmo shock RFP')
xlabel('time')
ylabel('Nuclear Localization')



%YFP Channel
figure(2)
    
fname_save = 'processed_data_Osmo_61_69.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW61_69 = all_tracks_vec;
all_times_vec_comb.EW61_69 = all_times_vec;

data_file_list = {'BMH37','EW61_69'};
clear('all_tracks_vec');
clear('all_times_vec');
all_tracks_vec = {};
all_times_vec = {};
for jj = 1:length(data_file_list);
    all_tracks_vec = [all_tracks_vec , all_tracks_vec_comb.(data_file_list{jj})];
    all_times_vec = [all_times_vec , all_times_vec_comb.(data_file_list{jj})];
end


clf 
hold on
channel = 'YFP'
legend_vec_YFP = {'BMH37 0.5M','BMH37 0.25M','EW61 0.5M', 'EW61 0.5M', 'EW69 0.5M', 'EW69 0.25M'}
cmap_YFP = [ 0,0,0;  %Black
    0.5,0.5,0.5;
    0,0,1;  %Blue
    0,1,1 %Blue 1
    0,1,0;  %Green
    0,0.75,0; %Green 1
           ];  
perm = [1,2,3,4,5,6]
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
title('Osmo shock YFP')
xlabel('time')
ylabel('Nuclear Localization')




%% Glucose Dropout plots
%RFP channel 


base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150729_EW049_61_69_check_GD\'

phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0,12]    


figure(3)
clf
hold on

channel = 'RFP'
legend_vec_RFP = {'BMH37 GD','BMH37 GD 0.1%','EW49 2 GD'}
cmap_RFP = [ 0,0,0;  %Black
    0.5,0.5,0.5; %Grey
    1,0,0  %Red
    ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3]  ;
N_RFP = length(perm)  
    
fname_save = 'processed_data_GD_37.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.BMH37 = all_tracks_vec;
all_times_vec_comb.BMH37 = all_times_vec;

fname_save = 'processed_data_GD_49.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW49 = all_tracks_vec;
all_times_vec_comb.EW49 = all_times_vec;

data_file_list = {'BMH37','EW49'};
clear('all_tracks_vec');
clear('all_times_vec');
all_tracks_vec = {};
all_times_vec = {};
for jj = 1:length(data_file_list);
    all_tracks_vec = [all_tracks_vec , all_tracks_vec_comb.(data_file_list{jj})];
    all_times_vec = [all_times_vec , all_times_vec_comb.(data_file_list{jj})];
end



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
title('Glucose Dropout RFP')
xlabel('time')
ylabel('Nuclear Localization')



%YFP Channel
figure(4)
    
fname_save = 'processed_data_GD_61_69.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW61_69 = all_tracks_vec;
all_times_vec_comb.EW61_69 = all_times_vec;

data_file_list = {'BMH37','EW61_69'};
clear('all_tracks_vec');
clear('all_times_vec');
all_tracks_vec = {};
all_times_vec = {};
for jj = 1:length(data_file_list);
    all_tracks_vec = [all_tracks_vec , all_tracks_vec_comb.(data_file_list{jj})];
    all_times_vec = [all_times_vec , all_times_vec_comb.(data_file_list{jj})];
end


clf 
hold on
channel = 'YFP'
legend_vec_YFP = {'BMH37 GD','BMH37 GD 0.1%','EW61 GD', 'EW61 GD 0.1%', 'EW69 GD', 'EW69 GD 0.1%'}
cmap_YFP = [ 0,0,0;  %Black
    0.5,0.5,0.5;
    0,0,1;  %Blue
    0,1,1 %Blue 1
    0,1,0;  %Green
    0,0.75,0; %Green 1
           ];  
perm = [1,2,3,4,5,6]

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
title('Glucose Dropout YFP')
xlabel('time')
ylabel('Nuclear Localization')

return

%% 1-NM-PP1 plots

%RFP channel 


base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150727_EW050_62_70_check_NMPP1\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    


figure(5)
clf
hold on

channel = 'RFP'
legend_vec_RFP = {'BMH42','EW50 1','EW50 2'}
cmap_RFP = [ 0,0,0;  %Black
    1,0,0;  %Red
    0.5,0,0 %Red 2
        ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3]  ;
N_RFP = length(perm)  
    
fname_save = 'processed_data_NMPP1_42.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.BMH42 = all_tracks_vec;
all_times_vec_comb.BMH42 = all_times_vec;

fname_save = 'processed_data_NMPP1_50.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW50 = all_tracks_vec;
all_times_vec_comb.EW50 = all_times_vec';

data_file_list = {'BMH42','EW50'};
clear('all_tracks_vec');
clear('all_times_vec');
all_tracks_vec = {};
all_times_vec = {};
for jj = 1:length(data_file_list);
    all_tracks_vec = [all_tracks_vec , all_tracks_vec_comb.(data_file_list{jj})];
    all_times_vec = [all_times_vec ; all_times_vec_comb.(data_file_list{jj})];
end



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
title('4uM NMPP1 RFP')
xlabel('time')
ylabel('Nuclear Localization')



%YFP Channel
figure(6)
    
fname_save = 'processed_data_NMPP1_62_70.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW62_70 = all_tracks_vec;
all_times_vec_comb.EW62_70 = all_times_vec';

data_file_list = {'BMH42','EW62_70'};
clear('all_tracks_vec');
clear('all_times_vec');
all_tracks_vec = {};
all_times_vec = {};
for jj = 1:length(data_file_list);
    all_tracks_vec = [all_tracks_vec , all_tracks_vec_comb.(data_file_list{jj})];
    all_times_vec = [all_times_vec ; all_times_vec_comb.(data_file_list{jj})];
end


clf 
hold on
channel = 'YFP'
legend_vec_YFP = {'BMH42','EW62 1','EW62 2', 'EW62 3', 'EW70 1', 'EW70 2', 'EW70 3'}
cmap_YFP = [ 0,0,0;  %Black
    0,0,1;  %Blue
    0,1,1 %Blue 1
    0,0.2,0.5 %Blue 2
    0,1,0;  %Green
    0,0.75,0; %Green 1
    0,0.5,0; %Green 2
        ];  
perm = [1,2,3,4,5,6,7]
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
title('4uM NMPP1 YFP')
xlabel('time')
ylabel('Nuclear Localization')




return

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


