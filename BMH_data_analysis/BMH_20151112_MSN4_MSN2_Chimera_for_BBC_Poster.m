function [] = BMH_20151112_MSN4_MSN2_Chimera_for_BBC_Poster()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)


%% 

%Osmo, NMPP1 for 
%Figure 1: SC.Msn2
%Figure 2: KL.Msn2
%Figure 3: SC.Msn2(KL.NLS)
%Figure 4: KL.Msn2(SC.NLS)
%Figure 5: SC.Msn4

%Osmo wells: {'D2','C2','F2','E2','B2'}
%NMPP1 wells: {'D10','C10','F10','E10','B10'}
%SDC wells: {'H10','G10','B11','A11','F10'}

%Load Osmo experiment
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20151014_47_52_54_60_76_Osmo_p5_p25\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    



% Sorb
% Row #	Strain	Media 		YFP channel*	
% B2 	BMH47y	0.5M		SC.MSN4         
% C2	EW52	"       	KL.MSN2         
% D2	EW54	"       	SC.MSN2         
% E2	EW60	"       	KL.MSN2(SC.NLS)	
% F2	EW76	"       	SC.MSN2(KL.NLS)	
% B3	BMH47y	0.25M   	SC.MSN4	
% C3	EW52	"       	KL.MSN2	
% D3	EW54	"       	SC.MSN2	
% E3	EW60	"       	KL.MSN2(SC.NLS)	
% F3	EW76	"       	SC.MSN2(KL.NLS)	

    
fname_save = 'processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_times_vec_osmo = all_times_vec;
all_tracks_vec_osmo = all_tracks_vec;

%Load NMPP1 experiment
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20151008_BMH47_EW52_54_60_76_NMPP1_SDC\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    


% NMPP1
% Row #	Strain	Media 		YFP channel*	
% A10 	BMH47y	4uM 1NMPP1	SC.MSN4         
% B10	EW52	"       	KL.MSN2         
% C10	EW54	"       	SC.MSN2         
% D10	EW60	"       	KL.MSN2(SC.NLS)	
% E10	EW76	"       	SC.MSN2(KL.NLS)	
% F10	BMH47y	SDC     	SC.MSN4	
% G10	EW52	"       	KL.MSN2	
% H10	EW54	"       	SC.MSN2	
% A10	EW60	"       	KL.MSN2(SC.NLS)	
% B10	EW76	"       	SC.MSN2(KL.NLS)	
    
fname_save = 'processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_times_vec_nmpp1 = all_times_vec;
all_tracks_vec_nmpp1 = all_tracks_vec;


%combine top osmo conditions and NMPP1 and SDC conditions

all_tracks_vec = [{all_tracks_vec_osmo{1:5}},all_tracks_vec_nmpp1]
all_times_vec = [{all_times_vec_osmo{1:5}},all_times_vec_nmpp1]

%% RFP control 

channel = 'RFP'
legend_vec_RFP = {'BMH47y 0.5M Sorb', 'EW52 0.5M Sorb', 'EW54 0.5M Sorb' ,'EW60 0.5M Sorb','EW76 0.5M Sorb','BMH47y 4uM 1NMPP1', 'EW52 4uM 1NMPP1', 'EW54 4uM 1NMPP1' ,'EW60 4uM 1NMPP1','EW76 4uM 1NMPP1','BMH47y SDC', 'EW52 SDC', 'EW54 SDC' ,'EW60 SDC','EW76 SDC'}
cmap_RFP = [  1,0,0;  %red
    1,0,0;  %red
    1,0,0;  %red
    1,0,0;  %red
    1,0,0;  %red
    0,0,0.5;  %light blue
    0,0,0.5;  %light blue
    0,0,0.5;  %light blue
    0,0,0.5;  %light blue
    0,0,0.5;  %light blue
    0,0,0;  %black
    0,0,0;  %black
    0,0,0;  %black
    0,0,0;  %black
    0,0,0  %black
          ];  

linestyle_vec_RFP = {'-','--','-.',':','-','-','--','-.',':','-','-','--','-.',':','-'}
marker_vec_RFP = {'None','None','None','None','x','None','None','None','None','x','None','None','None','None','x'}

%SC.MSN2 for all strains
figure(1)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [3,2,5,4,1,8,7,10,9,6,13,12,15,14,11];
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
axis([0,90,1,70])
title('SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')

%%YFP Channel

channel = 'YFP'
legend_vec_YFP = {'SC.MSN4 0.5M Sorb','KL.MSN2 0.5M Sorb','SC.MSN2 0.5M Sorb','KL.MSN2(SC.NLS) 0.5M Sorb','SC.MSN2(KL.NLS) 0.5M Sorb','SC.MSN4 4uM 1NMPP1','KL.MSN2 4uM 1NMPP1','SC.MSN2 4uM 1NMPP1','KL.MSN2(SC.NLS) 4uM 1NMPP1','SC.MSN2(KL.NLS) 4uM 1NMPP1', 'SC.MSN4 SDC','KL.MSN2 SDC','SC.MSN2 SDC','KL.MSN2(SC.NLS) SDC','SC.MSN2(KL.NLS) SDC'}
cmap_YFP = cmap_RFP

linestyle_vec_YFP = linestyle_vec_RFP
marker_vec_YFP = marker_vec_RFP

figure(2)
clf
hold on

perm = [3,2,5,4,1,8,7,10,9,6,13,12,15,14,11];
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



%Load NMPP1 experiment 





%% EXP 1 of 2 - EW47, 52, 54, 60 and 76 Osmo shock 0.5 and 0.25M
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20151014_47_52_54_60_76_Osmo_p5_p25\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    



% Sorb
% Row #	Strain	Media 		YFP channel*	
% B2 	BMH47y	0.5M		SC.MSN4         
% C2	EW52	"       	KL.MSN2         
% D2	EW54	"       	SC.MSN2         
% E2	EW60	"       	KL.MSN2(SC.NLS)	
% F2	EW76	"       	SC.MSN2(KL.NLS)	
% B3	BMH47y	0.25M   	SC.MSN4	
% C3	EW52	"       	KL.MSN2	
% D3	EW54	"       	SC.MSN2	
% E3	EW60	"       	KL.MSN2(SC.NLS)	
% F3	EW76	"       	SC.MSN2(KL.NLS)	

    
fname_save = 'processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')



%% RFP control 

channel = 'RFP'
legend_vec_RFP = {'BMH47y 0.25M Sorb', 'EW52 0.25M Sorb', 'EW54 0.25M Sorb' ,'EW60 0.25M Sorb','EW76 0.25M Sorb','BMH47y 0.5M Sorb', 'EW52 0.5M Sorb', 'EW54 0.5M Sorb' ,'EW60 0.5M Sorb','EW76 0.5M Sorb'}
cmap_RFP = [  0.5,0,0;  %light red
    0.5,0,0;  %light red
    0.5,0,0;  %light red
    0.5,0,0;  %light red
    0.5,0,0;  %light red
    1,0,0;  %red
    1,0,0;  %red
    1,0,0;  %red
    1,0,0;  %red
    1,0,0;  %red
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
legend_vec_YFP = {'SC.MSN4 0.25M Sorb','KL.MSN2 0.25M Sorb','SC.MSN2 0.25M Sorb','KL.MSN2(SC.NLS) 0.25M Sorb','SC.MSN2(KL.NLS) 0.25M Sorb', 'SC.MSN4 0.5M Sorb','KL.MSN2 0.5M Sorb','SC.MSN2 0.5M Sorb','KL.MSN2(SC.NLS) 0.5M Sorb','SC.MSN2(KL.NLS) 0.5M Sorb'}
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


%}

%{
%% EXP 2 of 2 - EW47, 52, 54, 60 and 76 Glucose Dropout full and 0.1%
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20151014_47_52_54_60_76_GD_full_p1\'

phases =  {'Pre', 'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0, 3.25+9]    



% GD
% Row #	Strain	Media 		YFP channel*	
% D6 	BMH47y	GD full		SC.MSN4         
% E6	EW52	"       	KL.MSN2         
% F6	EW54	"       	SC.MSN2         
% G6	EW60	"       	KL.MSN2(SC.NLS)	
% H6	EW76	"       	SC.MSN2(KL.NLS)	
% D7	BMH47y	GD 0.1%   	SC.MSN4	
% E7	EW52	"       	KL.MSN2	
% F7	EW54	"       	SC.MSN2	
% G7	EW60	"       	KL.MSN2(SC.NLS)	
% H7	EW76	"       	SC.MSN2(KL.NLS)	

    
fname_save = 'processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')



%% RFP control 

channel = 'RFP'
legend_vec_RFP = {'BMH47y GD full', 'EW52 GD full', 'EW54 GD full' ,'EW60 GD full','EW76 GD full','BMH47y GD 0.1%', 'EW52 GD 0.1%', 'EW54 GD 0.1%' ,'EW60 GD 0.1%','EW76 GD 0.1%'}
cmap_RFP = [  0,0,1;  %Blue
     0,0,1;  %Blue
     0,0,1;  %Blue
     0,0,1;  %Blue
     0,0,1;  %Blue
     0,0,0.5;  %light Blue
     0,0,0.5;  %light Blue
     0,0,0.5;  %light Blue
     0,0,0.5;  %light Blue
     0,0,0.5;  %light Blue
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
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
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
legend_vec_YFP = {'SC.MSN4 GD full','KL.MSN2 GD full','SC.MSN2 GD full','KL.MSN2(SC.NLS) GD full','SC.MSN2(KL.NLS) GD full', 'SC.MSN4 GD 0.1%','KL.MSN2 GD 0.1%','SC.MSN2 GD 0.1%','KL.MSN2(SC.NLS) GD 0.1%','SC.MSN2(KL.NLS) GD 0.1%'}
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
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_YFP_plot) %,'Location','NE');
axis([0,70,1,8])
title('YFP Channel')
xlabel('time')
ylabel('Nuclear Localization')

%}




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


