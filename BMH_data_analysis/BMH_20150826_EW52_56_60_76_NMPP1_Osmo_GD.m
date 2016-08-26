function [] = BMH_20150826_EW52_56_60_76_NMPP1_Osmo_GD()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%{
%% EXP 1 of 2 - EW52, 56, 60 and 76 4uM NMPP1 0.5M Osmo shock
base_dir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\20150826_EW52_56_60_76_Osmo_NMPP1_GD\NMPP1_Osmo\'

phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    


%Exp 1
%all locations had 4 sites
% Osmo shock 0.5M Sorbitol added after Pre phase (approx 10 min)
% Well# Strain	fluor	Condition
% A7	52      both	4uM NMPP1
% B7	56      YFP     4uM NMPP1
% C7	60      both	4uM NMPP1
% D7	76      both	4uM NMPP1
% E7	52      both	0.5M Sorb
% F7	56      YFP     0.5M Sorb
% G7	60      both	0.5M Sorb
% H7	76      both	0.5M Sorb

    
fname_save = 'processed_data_52_60_76.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')



%% RFP control 

channel = 'RFP'
legend_vec_RFP = {'EW52 4uM NMPP1','EW60 4uM NMPP1','EW76 4uM NMPP1', 'EW52 0.5M Sorb','EW60 0.5M Sorb','EW76 0.5M Sorb'}
cmap_RFP = [ 0,0.5,1;  %light blue
    0,.5,1; %light blue
    0,.5,1; %light blue
    1,0,0;  %Red
    1,0,0;  %Red
    1,0,0;  %Red
      ];  

linestyle_vec_RFP = {'-','--','-.','-','--','-.'}

%SC.MSN2 for all strains
figure(1)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3,4,5,6];
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
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
axis([0,70,1,7])
title('SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')


%%YFP Channel

fname_save = 'processed_data_52_60_76.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW52_60_76 = all_tracks_vec;
all_times_vec_comb.EW52_60_76 = all_times_vec;

fname_save = 'processed_data_56.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW56 = all_tracks_vec;
all_times_vec_comb.EW56 = all_times_vec;

data_file_list = {'EW52_60_76','EW56'};
clear('all_tracks_vec');
clear('all_times_vec');
all_tracks_vec = {};
all_times_vec = {};
for jj = 1:length(data_file_list);
    all_tracks_vec = [all_tracks_vec , all_tracks_vec_comb.(data_file_list{jj})];
    all_times_vec = [all_times_vec , all_times_vec_comb.(data_file_list{jj})];
end



channel = 'YFP'
legend_vec_YFP = {'KL.MSN2 4uM NMPP1','KL.MSN2(SC.NLS) 4uM NMPP1','SC.MSN2(KL.NLS) 4uM NMPP1', 'KL.MSN2 0.5M Sorb','KL.MSN2(SC.NLS) 0.5M Sorb','SC.MSN2(KL.NLS) 0.5M Sorb', 'SC.MSN2 4uM NMPP1', 'SC.MSN2 0.5M Sorb' }
cmap_YFP = [ 0,0.5,1;  %light blue
    0,.5,1; %light blue
    0,.5,1; %light blue
    1,0,0;  %Red
    1,0,0;  %Red
    1,0,0;  %Red
    0,.5,1; %light blue
    1,0,0;  %Red
    ];  

linestyle_vec_YFP = {'-','--','-.','-','--','-.',':',':'}

%EW52 (KL.MSN2) all perturbations
figure(2)
clf
hold on

perm = [1,4];
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
figure(3)
clf
hold on

perm = [2,5];
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

perm = [3,6];
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

perm = [7,8];
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
title('SC.MSN2-YFP')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')


figure(6)
clf
hold on

perm = [1,2,3,7];
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
title('NMPP1')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')

figure(7)
clf
hold on

perm = [4,5,6,8];
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
title('Osmo shock')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')

%}


%% EXP 2 of 2 - EW52, 56, 60 and 76 Glucose Dropout
base_dir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\20150826_EW52_56_60_76_Osmo_NMPP1_GD\GD\'

phases =  {'Pre', 'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0,10.75]    


%Exp 1
%all locations had 4 sites
% Osmo shock 0.5M Sorbitol added after Pre phase (approx 10 min)
% Well# Strain	fluor	Condition
% A12	52      both	GD 
% B12   56      YFP     GD 
% C12	60      both	GD
% D12	76      both	GD
% E12	52      both	GD 0.1%
% F12	56      YFP     GD 0.1%
% G12	60      both	GD 0.1%
% H12	76      both	GD 0.1%

    
fname_save = 'processed_data_52_60_76.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')


%% RFP control 

channel = 'RFP'
legend_vec_RFP = {'EW52 GD 0%','EW60 GD 0%','EW76 GD 0%', 'EW52 GD 0.1%','EW60 0.1%','EW76 0.1%'}
cmap_RFP = [ 0,0,1;  %blue
    0,0,1; %blue
    0,0,1; %blue
    0,0.75,1;  %light blue
    0,.75,1; %light blue
    0,.75,1; %light blue
      ];  

linestyle_vec_RFP = {'-','--','-.','-','--','-.'}

%SC.MSN2 for all strains
figure(8)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3,4,5,6];
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
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_RFP_plot) %,'Location','NE');
axis([0,70,1,7])
title('SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')


%%YFP Channel

fname_save = 'processed_data_52_60_76.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW52_60_76 = all_tracks_vec;
all_times_vec_comb.EW52_60_76 = all_times_vec;

fname_save = 'processed_data_56.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
all_tracks_vec_comb.EW56 = all_tracks_vec;
all_times_vec_comb.EW56 = all_times_vec;

data_file_list = {'EW52_60_76','EW56'};
clear('all_tracks_vec');
clear('all_times_vec');
all_tracks_vec = {};
all_times_vec = {};
for jj = 1:length(data_file_list);
    all_tracks_vec = [all_tracks_vec , all_tracks_vec_comb.(data_file_list{jj})];
    all_times_vec = [all_times_vec , all_times_vec_comb.(data_file_list{jj})];
end



channel = 'YFP'
legend_vec_YFP = {'KL.MSN2 GD 0%','KL.MSN2(SC.NLS) GD 0%','SC.MSN2(KL.NLS) GD 0%', 'KL.MSN2 GD 0.1%','KL.MSN2(SC.NLS) GD 0.1%','SC.MSN2(KL.NLS) GD 0.1%', 'SC.MSN2 GD 0%', 'SC.MSN2 GD 0.1%' }


cmap_YFP = [ 0,0,1;  %blue
    0,0,1; %blue
    0,0,1; %blue
    0,0.75,1;  %light blue
    0,0.75,1;  %light blue
    0,0.75,1;  %light blue
    0,0,1; %blue
    0,0.75,1;  %light blue
    ];  

linestyle_vec_YFP = {'-','--','-.','-','--','-.',':',':'}

%EW52 (KL.MSN2) all perturbations
figure(9)
clf
hold on

perm = [1,4];
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
figure(10)
clf
hold on

perm = [2,5];
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


figure(11)
clf
hold on

perm = [3,6];
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


figure(12)
clf
hold on

perm = [7,8];
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
title('SC.MSN2-YFP')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')


figure(13)
clf
hold on

perm = [1,2,3,7];
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
title('GD 0%')
axis([0,70,1,7])
xlabel('time')
ylabel('Nuclear Localization')

figure(14)
clf
hold on

perm = [4,5,6,8];
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
title('GD 0.1%')
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


