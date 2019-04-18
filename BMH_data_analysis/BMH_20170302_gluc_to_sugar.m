function [] = BMH_20170302_gluc_to_sugar()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)




%% NMPP1 Does Response
%base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20170302_KL_sugar_replacement_3\'
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20170302_SC_sugar_replacement_3\'

phases =  {'Pre', 'Post'}%,'Post_p2'} %,'Post'} 
%shift_timing = [0,10]    
shift_timing = [0,14]    

% Column 10 for KLac, 
% Column 10 for SC yEW051
% Column 11 for SC ySV076
% A12 is ySV076 + 0.5M Xylose osmo shock
%
% Well	condition	Results (visual)
% A: GD: 
% B: SDC: 
% C: Sorb: 
% D; Gal
% E: EtOH 2%
% F: EtOH 2/3%
% G: Glycerol 1.6%
% H: Glycerol 0.8%




fname_save = 'processed_data.mat';
%fname_save = 'processed_data_SV076.mat';

load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')


%% RFP channel 
%channels = {'RFP_wheel'};
channels = {'RFP_wheel','YFP_wheel'};
titles = {'SC.MSN2', 'KL.MSN2 in SC cells'};
%titles = {'KL.MSN2', 'SC.MSN2 in KL cells'};
%titles = {'Hog1'};

for kk = 1:length(channels)
    channel = channels{kk};
    legend_vec_RFP = {'H20','SDC', 'Sorb 0.11M', 'Gal 0.11M', '2% EtOH', '0.67% EtOH','1.6% Gly', '0.8% Gly'};
    %legend_vec_RFP = {'H20','SDC', 'Sorb 0.11M', 'Gal 0.11M', '2% EtOH', '0.67% EtOH','1.6% Gly', '0.8% Gly','2% glu + 0.5M Xylose'};
    cmap_RFP = [ 0, 0, 0;
    0, 0, 0;
    66,128,244; %blue
    80,244,66; %green
    244,66,238; %purple
    244,66,238; %purple
    244,122,66; %orange
    244,122,66 %orange
    %255,0,0  %red - only for Hog
    ]/255;

    %linestyle_vec_RFP = {'-','--','-','-','-','--','-','--','-'}
    %marker_vec_RFP = {'None','None','None','None','None','None','None','None','None'}

    linestyle_vec_RFP = {'-','--','-','-','-','--','-','--'}
    marker_vec_RFP = {'None','None','None','None','None','None','None','None'}

    
    %All Strains and doses
    figure(kk)
    clf
    hold on

    %cmap = jet(length(legend_vec_RFP));
    %perm = 1:length(legend_vec_RFP);
    %perm = [1,2,3,4,5,6,7,8,9];
    perm = [1,2,3,4,5,6,7,8];
    N_RFP = length(perm);

    legend_vec_RFP_plot = legend_vec_RFP(perm);
    cmap = cmap_RFP(perm,:);
    plt_grp = zeros(length(legend_vec_RFP_plot),1);
    for jj = 1:N_RFP
        all_tracks = all_tracks_vec{perm(jj)};
        all_times = all_times_vec{perm(jj)};
        color_val = cmap(jj,:);
        plt_grp(jj) = hggroup;
        plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_RFP{perm(jj)},'Marker',marker_vec_RFP{perm(jj)}};
        for ph = 1: length(phases)
            tracks = all_tracks.(phases{ph});
            timevals = all_times.(phases{ph});
            p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
            set(p,'Parent',plt_grp(jj))
        end
    end

    hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
    axis([0,60,1,8])
    title(titles{kk})
    xlabel('time')
    ylabel('Nuclear Localization')
end

return

%%YFP Channel

channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP

linestyle_vec_YFP = linestyle_vec_RFP
marker_vec_YFP = marker_vec_RFP

%SC.MSN2 for all strains
figure(2)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [2,4,6,1,3,5];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)},'Marker',marker_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
axis([0,130,1,3])
title('SC.MSN2-YFP in KL Cells')
xlabel('time')
ylabel('Nuclear Localization')



fname_save = 'processed_data_SC.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec;
all_times_vec;

%% RFP channel 

channel = 'RFP'
legend_vec_RFP = {'PKA AS, 4uM 1NMPP1 + 0.5M Sorb','PKA AS, 0.5M Sorb'}
cmap_RFP = [       1,0,1;  %Purple
     1,0,0;  %Red
          ];  

linestyle_vec_RFP = {'-','-'}
marker_vec_RFP = {'None','None'}

%KL.MSN2 for all strains
figure(3)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [2,1];
N_RFP = length(perm);

legend_vec_RFP_plot = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plt_grp = zeros(length(legend_vec_RFP_plot),1);
for jj = 1:N_RFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_RFP{perm(jj)},'Marker',marker_vec_RFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
axis([0,130,1,7])
title('SC.MSN2-RFP in SC Cells')
xlabel('time')
ylabel('Nuclear Localization')



%%YFP Channel

channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP

linestyle_vec_YFP = linestyle_vec_RFP
marker_vec_YFP = marker_vec_RFP

%SC.MSN2 for all strains
figure(4)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [2,1];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)},'Marker',marker_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
axis([0,130,1,7])
title('KL.MSN2-YFP in SC Cells')
xlabel('time')
ylabel('Nuclear Localization')





return



