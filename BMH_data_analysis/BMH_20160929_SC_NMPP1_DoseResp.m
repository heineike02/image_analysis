function [] = BMH_20160929_SC_NMPP1_DoseResp()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)




%% SC NMPP1 Does Response
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20160929_SC_NMPP1_Dose_Resp\'

phases =  {'Exp'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    


% 
% Well	Strain	Concentration	Results (visual)
% A1	AS	0	No loc
% B1	WT	0	
% C1	AS	10nM	
% D1	AS	25nM	
% E1	AS	100nM	
% F1	WT	100nM	
% A2	AS	250nM	Some Loc
% B2	AS	1uM	
% C2	AS	4uM	
% D2	WT 	4uM	
% E2	AS	10uM	
% F2	WT	10uM	
% G2	WT	0 10 DMSO: 490 SDC	


fname_save = 'processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec;
all_times_vec;

%% RFP channel 

channel = 'YFP'
legend_vec_RFP = {'AS Cont','WT 0', 'AS 10nM', 'AS 25 nM', 'AS 100nM', 'WT 100nM', 'AS 250nM', 'AS 1uM', 'AS 4uM', 'WT 4uM', 'AS 10uM', 'WT 10uM','WT DMSO'}
cmap_RFP = [ 0, 0, 0;
0, 0, 0;
67, 42, 27;
156, 97, 62;
253, 158, 100;
253, 158, 100;
254, 160, 102;
255, 161, 103;
255, 182, 116;
255, 182, 116;
255, 199, 127;
255, 199, 127;
0, 0, 0
]/255;

linestyle_vec_RFP = {'-','--','-','-','-','--','-','-','-','--','-','--','-.'}
marker_vec_RFP = {'None','None','None','None','None','None','None','None','None','None','None','None','None' }

%All Strains and doses
figure(1)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,3,4,5,7,8,9,11,2,13,6,10,12];
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
axis([0,25,1,8])
title('Dose Respones to NMPP1')
xlabel('time')
ylabel('Nuclear Localization')


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



