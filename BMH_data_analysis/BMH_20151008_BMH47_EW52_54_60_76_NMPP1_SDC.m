function [] = BMH_20151008_BMH47_EW52_54_60_76_NMPP1_SDC()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)


%% EXP 1 of 2 - EW52, 60 and 76 4uM NMPP1 0.5M Osmo shock
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20151008_BMH47_EW52_54_60_76_NMPP1_SDC\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    




% Row #	Strain	Media 	Vol	YFP channel*	Notes
% A10 	BMH47y	NMPP1	150	SC.MSN4         Msn4 localizes, but transiently returns from nucleus
% B10	EW52	"       "	KL.MSN2         Many cells washed away.  Some cells very nuclear, most only transient.  Worried there may be some contamination. 
% C10	EW54	"       "	SC.MSN2         Both localize intensely
% D10	EW60	"       "	KL.MSN2(SC.NLS)	
% E10	EW76	"       "	SC.MSN2(KL.NLS)	
% F10	BMH47y	SDC     "	SC.MSN4	
% G10	EW52	"       "	KL.MSN2	
% H10	EW54	"       "	SC.MSN2	
% A11	EW60	"       "	KL.MSN2(SC.NLS)	
% B11	EW76	"       "	SC.MSN2(KL.NLS)	

    
fname_save = 'processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')



%% RFP control 

channel = 'RFP'
legend_vec_RFP = {'BMH47y 4uM 1NMPP1', 'EW52 4uM 1NMPP1', 'EW54 4uM 1NMPP1' ,'EW60 4uM 1NMPP1','EW76 4uM 1NMPP1','BMH47y 4uM SDC', 'EW52 4uM SDC', 'EW54 4uM SDC' ,'EW60 4uM SDC','EW76 4uM SDC'}
cmap_RFP = [ 0,0.5,1;  %light blue
    0,0.5,1; %light blue
    0,0.5,1; %light blue
    0,0.5,1;  %light blue
    0,0.5,1;  %light blue
    0,0,0;  %Black
    0,0,0;  %Black
    0,0,0; %Black
    0,0,0; %Black
    0,0,0; %Black
       ];  

linestyle_vec_RFP = {'-','--','-.',':','-','-','--','-.',':','-'}
marker_vec_RFP = {'None','None','None','None','x','None','None','None','None','x'}

%SC.MSN2 for all strains
figure(1)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3,4,5,6,7,8,9,10];
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
axis([0,70,1,8])
title('SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')


%%YFP Channel


channel = 'YFP'
legend_vec_YFP = {'SC.MSN4 4uM 1NMPP1','KL.MSN2 4uM 1NMPP1','SC.MSN2 4uM 1NMPP1','KL.MSN2(SC.NLS) 4uM 1NMPP1','SC.MSN2(KL.NLS) 4uM 1NMPP1', 'SC.MSN4 4uM SDC','KL.MSN2 4uM SDC','SC.MSN2 4uM SDC','KL.MSN2(SC.NLS) 4uM SDC','SC.MSN2(KL.NLS) 4uM SDC'}
cmap_YFP = cmap_RFP

linestyle_vec_YFP = {'-','--','-.',':','-','-','--','-.',':','-'}
marker_vec_YFP = {'None','None','None','None','x','None','None','None','None','x'}

figure(2)
clf
hold on

perm = [1,2,3,4,5,6,7,8,9,10];
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
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
axis([0,70,1,8])
title('YFP Channel')
xlabel('time')
ylabel('Nuclear Localization')


return



%EW60 (KL.MSN2(SC.NLS)) all perturbations
figure(3)
clf
hold on

perm = [2,5,8];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
title('KL.MSN2(SC.NLS)-YFP')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')


figure(4)
clf
hold on

perm = [3,6,9];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
title('SC.MSN2(KL.NLS)-YFP')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')


figure(5)
clf
hold on

perm = [1,2,3,7,8,9];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
title('NMPP1 v.s. Control')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')

figure(6)
clf
hold on

perm = [4,5,6,7,8,9];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
title('Osmo shock v.s. Control')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')




%% EXP 2 of 2 - EW52, 60 and 76 Glucose Dropout
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150821_52_60_76_GD\'

phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0,11.75]    

%Strain BMH42 - SC.MSN2-RFP,  KL.MSN2-YFP 
%

%Exp 1
%all locations had 4 sites
% Glucose drop added after 9 min (3 images).  Perturbation took 2m45s
% Well#	Strain	fluor	Condition
% A11	52	both	GD
% B11	60	both	GD
% C11	76	both	GD
% D11	52	both	GD 0.1%
% E11	60	both	GD 0.1%
% F11	76	both	GD 0.1%
% G11	52	both	SDC
% H11	60	both	SDC
% C10	76	Both	SDC

    
fname_save = 'processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')


%% RFP control 

channel = 'RFP'
legend_vec_RFP = {'EW52 GD 0%','EW60 GD0%','EW76 GD0%', 'EW52 GD 0.1%','EW60 GD 0.1%','EW76 0.2%', 'EW52 SDC','EW60 SDC','EW76 SDC' }
cmap_RFP = [ 0,0,1;  %blue
    0,0,1; %blue
    0,0,1; %blue
    0,0.75,1;  %light blue
    0,.75,1; %light blue
    0,.75,1; %light blue
    0,0,0; %Black
    0,0,0; %Black
    0,0,0; %Black
       ];  

linestyle_vec_RFP = {'-','--','-.','-','--','-.','-','--','-.'}

%SC.MSN2 for all strains
figure(7)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3,4,5,6,7,8,9];
N_RFP = length(perm);

legend_vec_RFP_plot = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plt_grp = zeros(length(legend_vec_RFP_plot),1);
for jj = 1:N_RFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_RFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
axis([0,70,1,7])
title('SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')


%%YFP Channel


channel = 'YFP'
legend_vec_YFP = {'KL.MSN2 GD 0%','KL.MSN2(SC.NLS) GD 0%','SC.MSN2(KL.NLS) GD 0%', 'KL.MSN2 GD 0.1%','KL.MSN2(SC.NLS) GD 0.1%','SC.MSN2(KL.NLS) GD 0.1%', 'KL.MSN2 SDC','KL.MSN2(SC.NLS) SDC','SC.MSN2(KL.NLS) SDC' }
cmap_YFP = cmap_RFP

linestyle_vec_YFP = {'-','--','-.','-','--','-.','-','--','-.'}

%EW52 (KL.MSN2) all perturbations
figure(8)
clf
hold on

perm = [1,4,7];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
axis([0,70,1,7])
title('KL.MSN2-YFP')
xlabel('time')
ylabel('Nuclear Localization')


%EW60 (KL.MSN2(SC.NLS)) all perturbations
figure(9)
clf
hold on

perm = [2,5,8];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
title('KL.MSN2(SC.NLS)-YFP')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')


figure(10)
clf
hold on

perm = [3,6,9];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
title('SC.MSN2(KL.NLS)-YFP')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')


figure(11)
clf
hold on

perm = [1,2,3,4,5,6,7,8,9];
N_YFP = length(perm);

legend_vec_YFP_plot = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plt_grp = zeros(length(legend_vec_YFP_plot),1);
for jj = 1:N_YFP
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_YFP{perm(jj)}};
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
title('Glucose Dropout')
axis([0,70,1,7])
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


