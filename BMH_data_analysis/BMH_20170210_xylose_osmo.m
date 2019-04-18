function [] = BMH_20170210_xylose_osmo()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)




%% NMPP1 Does Response
%base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20170210_KL_xylose_osmo\'
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20170210_SC_xylose_osmo\'

phases =  {'Pre', 'Post'}%,'Post_p2'} %,'Post'} 
%shift_timing = [0,10]  %KL
shift_timing = [0,14]  %SC  

% Column 9 for KLac, 
% Column 7 for SC yEW051
% Column 8 for SC ySV076
%
%
% Well	condition	Results (visual)
% A: 2% Gluc + 0.5M Xyl 
% B: 2% Gluc + 0.5M Sorb
% C: 0.11M Xyl
% D; 0.11M Sorb
% E: H20
% F: 0.11M Xyl + 0.5M Xyl
% G: 0.11M Sorb + 0.5M Sorb
% H: 2% Gluc

%fname_save = 'processed_data.mat';
%fname_save = 'processed_data_Hog1.mat';
fname_save = 'processed_data_Msn2.mat';
[base_dir,fname_save]
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')


%% RFP channel 
channels = {'RFP_wheel','YFP_wheel'};
%channels = {'RFP_wheel'}
%titles = {'KL.MSN2', 'SC.MSN2 in KL cells'};
titles = {'SC.MSN2', 'KL.MSN2 in SC cells'};
%titles = {'Hog1'};

for kk = 1:length(channels)
    channel = channels{kk};
    legend_vec_RFP = {'Osmo Shock: 2% Gluc + 0.5M Xyl', '2% Gluc + 0.5M Sorb', 'Glucose Dropout: Osmo balance w/0.11M Xyl','0.11M Sorb','Glucose Dropuot: no Osmo balance','0.11M Xyl + 0.5M Xyl','0.11M Sorb + 0.5M Sorb','2% Gluc'}
    cmap_RFP = [255,0,0; %Red
    66,128,244; %blue
    255,0,0; %Red
    66,128,244; %blue
    0,0,0; 
    255,0,0; %Red
    66,128,244; %blue
    0,0,0
    ]/255;

    linestyle_vec_RFP = {'--','--','-','-','-',':',':','--'}
    marker_vec_RFP = {'None','None','None','None','None','None','None','None'}

    %All Strains and doses
    figure(kk)
    clf
    hold on

    %cmap = jet(length(legend_vec_RFP));
    %perm = 1:length(legend_vec_RFP);
    %perm = [1,2,3,4,5,6,7,8];
    perm = [1,3,5,8];
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
    axis([0,60,1,8]) %SC Msn2
    %axis([0,60,1,3]) %SC Hog1
    %axis([0,90,1,4]) %KL Msn2
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



