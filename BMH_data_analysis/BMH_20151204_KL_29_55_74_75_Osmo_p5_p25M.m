function [] = BMH_20151215_KL_NMPP1_Osmo()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)




%% 0.5M and 0.25M Osmo shock, K.Lactis
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20151103_29_55_74_75_osmo\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    



% Well	Strain	Condition	RFP	YFP
% A5	29	0.5M Osmo	KL.MSN2	NA
% B5	55 (29 -URA)	"	KL.MSN2	NA
% C5	74	"	KL.MSN2	SC.MSN2
% D5	75 (TPK2 AS)	"	KL.MSN2	SC.MSN2
% E5	29	0.25M Osmo	KL.MSN2	NA
% F5	55	"	KL.MSN2	NA
% G5	74	"	KL.MSN2	SC.MSN2
% H5	75	"	KL.MSN2	SC.MSN2


fname_save = 'processed_data_74_75_RFP_YFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec_2Color = all_tracks_vec;
all_times_vec_2Color = all_times_vec;

fname_save = 'processed_data_29_55_RFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec_RFP = all_tracks_vec;
all_times_vec_RFP = all_times_vec;

all_tracks_vec = [all_tracks_vec_2Color, all_tracks_vec_RFP];
all_times_vec = [all_times_vec_2Color, all_times_vec_RFP];

%% RFP channel 

channel = 'RFP'
legend_vec_RFP = {'74 0.5M Sorb','75 (TPK2-AS) 0.5M Sorb' ,'74 0.25M Sorb','75 (TPK2-AS) 0.25M Sorb', '29 0.5M Sorb', '55 (-URA) 0.5M Sorb','29 0.25M Sorb', '55 (-URA) 0.25M Sorb'}
cmap_RFP = [       1,0,0;  %Red
     1,0,0;  %Red
     0.5,0,0;  %Light Red
     0.5,0,0;  %Light Red
     1,0,0;  %Red
     1,0,0;  %Red
     0.5,0,0;  %Light Red
     0.5,0,0  %Light Red
          ];  

linestyle_vec_RFP = {'-','-.','-','-.','--',':','--',':'}
marker_vec_RFP = {'None','None','None','None','None','None','None','None'}

%KL.MSN2 for all strains
figure(1)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [5,6,1,2,7,8,3,4];
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
axis([0,110,1,3])
title('KL.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')



%%YFP Channel

fname_save = 'processed_data_74_75_RFP_YFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

channel = 'YFP'
legend_vec_YFP =  {'74 0.5M Sorb','75 (TPK2-AS) 0.5M Sorb', '74 0.25M Sorb','75 (TPK2-AS) 0.25M Sorb'}

cmap_YFP = [       1,0,0;  %Red
     1,0,0;  %Red
     0.5,0,0;  %Light Red
     0.5,0,0  %Light Red
          ];  
      
linestyle_vec_YFP = {'--',':','--',':'}
marker_vec_YFP = {'None','None','None','None',}

figure(2)
clf
hold on

perm = [1,2,3,4];
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
axis([0,110,1,3])
title('SC.Msn2-YFP')
xlabel('time')
ylabel('Nuclear Localization')


return



