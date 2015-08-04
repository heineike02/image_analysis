function [] = BMH_20150727_EW50_62_70_checks()
%

profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150727_EW050_62_70_check_osmo\'

phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0,11.5]    


%Strain BMH42 - SC.MSN2-RFP,  KL.MSN2-YFP 
%

%Exp 1
%all locations had 4 sites
% Osmo shock 0.5M Sorbitol added after Pre phase (approx 10 min)
%
% Glucose Dropout (0.1 % from 2%) followed by Osmo shock (0.35M) at various times.
% osmo balanced with sorbitol relative to osmotic stress condition.
%
%A1	 BMH42
%B1	 EW50 1
%C1	 EW50 2
%D1	 EW62 1
%E1	 EW62 2
%F1	 EW62 3
%G1	 EW70 1
%H1	 EW70 2
%A2  EW70 3 




%% Osmo shock plots

%RFP channel 

figure(1)
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
    
fname_save = 'processed_data_42.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.BMH42 = all_tracks_vec;
all_times_vec_comb.BMH42 = all_times_vec;

fname_save = 'processed_data_50.mat';
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
title('Osmo shock RFP')
xlabel('time')
ylabel('Nuclear Localization')



%YFP Channel
figure(2)
    
fname_save = 'processed_data_62_70.mat';
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
title('Osmo shock YFP')
xlabel('time')
ylabel('Nuclear Localization')



%Everything below return is not plotted
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


