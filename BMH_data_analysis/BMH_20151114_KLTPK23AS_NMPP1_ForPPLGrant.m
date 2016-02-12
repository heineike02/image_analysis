function [] = BMH_20151114_KLTPK23AS_NMPP1_ForPPLGrant()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)




%% NMPP1 for BMH 55, 75, 77 and 78  
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20151113_55_75_78prepost_nmpp1\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    



% Sorb
% RFP Channel is KL.MSN2
%
% Row #	Strain	Media 		YFP channel   AS mutants	Crispr plas
% A1 	BMH55y	4uM 1NMPP1	NA              NA          -
% B1	BMH75y	"       	SC.MSN2         TPK2        -
% C1	BMH77y	"       	"               TPK2/3      +
% D1	BMH78y	"       	"               "           -
% E1	BMH55y	SDC       	NA              NA          -
% F1	BMH75y	"          	SC.MSN2         TPK2        -
% G1	BMH77y	"       	"               TPK2/3      +
% H1	BMH78y	"       	"               TPK2/3      -


fname_save = 'processed_data_75_78.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec_75_78 = all_tracks_vec;
all_times_vec_75_78 = all_times_vec;

fname_save = 'processed_data_55.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec_55 = all_tracks_vec;
all_times_vec_55 = all_times_vec;

all_tracks_vec = [all_tracks_vec_75_78, all_tracks_vec_55];
all_times_vec = [all_times_vec_75_78, all_times_vec_55];

%% RFP channel 

channel = 'RFP'
legend_vec_RFP = {'TPK2-AS', 'TPK2/3-AS +Cas9','TPK2/3-AS','TPK2-AS SDC','TPK2/3-AS +Cas9 SDC','TPK2/3-AS SDC','WT', 'WT SDC'}
cmap_RFP = [       0,0,0.5;  %light Blue
     0,0,0.5;  %light Blue
     0,0,0.5;  %light Blue
    0,0,0;  %black
    0,0,0;  %black
    0,0,0;  %black
    0,0,0.5;  %light Blue
    0,0,0;  %black
          ];  

linestyle_vec_RFP = {'--','-.','-','--','-.','-',':',':'}
marker_vec_RFP = {'None','None','None','None','None','None','None','None'}

%KL.MSN2 for all strains
figure(1)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [3,7,6,8];
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
axis([0,100,1,3])
title('KL.MSN2-RFP Nuclear Localization, 4uM NMPP1')
xlabel('time')
ylabel('Nuclear Localization')


return


%KL.MSN2 for 55, 75, and 78 strains
figure(2)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [7,3,1,8,5,4];
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
axis([0,100,1,3])
title('KL.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')





%%YFP Channel

fname_save = 'processed_data_75_78.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

channel = 'YFP'
legend_vec_YFP = {'TPK2-AS 4uM NMPP1', 'TPK2/3-AS +Cas9 4uM NMPP1','TPK2/3-AS 4uM NMPP1','TPK2-AS SDC','TPK2/3-AS +Cas9 SDC','TPK2/3-AS SDC'}

cmap_YFP = [       0,0,0.5;  %light Blue
     0,0,0.5;  %light Blue
     0,0,0.5;  %light Blue
    0,0,0;  %black
    0,0,0;  %black
    0,0,0  %black
          ]; 
      
linestyle_vec_YFP = {'--','-.','-','--','-.','-'}
marker_vec_YFP = {'None','None','None','None','None','None'}

figure(3)
clf
hold on

perm = [2,3,1,5,6,4];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}, 'Marker',marker_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
axis([0,100,1,5])
title('SC.Msn2-YFP')
xlabel('time')
ylabel('Nuclear Localization')


figure(4)
clf
hold on

perm = [3,1,6,4];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}, 'Marker',marker_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
axis([0,100,1,5])
title('SC.Msn2-YFP')
xlabel('time')
ylabel('Nuclear Localization')


return



