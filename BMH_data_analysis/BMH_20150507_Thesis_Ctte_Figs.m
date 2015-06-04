function all_tracks = BMH_20150507_Thesis_Ctte_Figs()

ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

%Plots for thesis Ctte Mtg 06May15


%% Figs 1 and 2 - GD and Osmo Shock in 2-color exp
%Data from 20140716_39Y_sorb_GD_2Color
%all locations had 4 sites
%A8 SC 37_1 GD 
%B8 SC 37_2 GD 
%C8 SC 38_1 GD 
%D8 SC 38_2 GD 
%E8 SC 11-38 GD 
%F8 KL 39 GD  
%G8 SC 37_1 SDC to SDC+0.5M Sorb


phases =  {'Pre','Post'} %,'Post'} 
shift_timing = [0,8.66666]  
base_dir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\20140716\'

figure(1)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'Gluc Drop, 0.111M sorb';
    'strain 2 Gluc Drop, 0.111M sorb';
    '0.5M sorb'}
    
fname_save = '20140716_processed_data_SC_2color.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [0,0,1;  %Blue
    0,0.5,1;  %Light blue
    1,0,0;  %Red
    ];

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,3];
N_RFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,70,1,8])



figure(2)
clf
hold on
channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP;

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,3];
N_YFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_YFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, KL.MSN2-YFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,70,1,4])


%% Figs 3 and 4 1-NM-PP1 2-color

%data from 20140918_AS_2Color
%all locations had 4 sites
%4uM 1-NM-PP1 
%A1 37, 
%B1 42-1, 
%C1 42-2, 
%0.2uM: 
%D1 42-1, 
%E1 42-2 
%0 
%F1 37 
%G1 42-1, 
%H1 42-2

%Remove various sites from list 
%C1 site 3 starts off with no cells - Remove
%B1 site 3 starts off with no cells - Remove
%A1 Site 3 starts off with one cell
%posvecSP.SC{2,4} = 'NA'
%posvecSP.SC{3,4} = 'NA'

legend_vec = {'37 4uM','42_1 4uM', '42_2 4uM','42_1 0.2 uM1' , '42_2 0.2 uM', '37 cont','42_1 cont', '42_2 cont'} %'39y 1 0.5M sorb','39y 1 No gluc, 0.111M sorb','39y 2 0.5M sorb','39y 2 No gluc, 0.111M sorb', '29y 1 0.5M sorb','29y 1 No gluc, 0.111M sorb'}

wellvecSP.SC = {'A1','B1','C1','D1', 'E1', 'F1', 'G1', 'H1'};  %RFP_only 'C8','D8','E8'
Nsites = 4;


phases =  {'Pre','Post'} %,'Post'} 
shift_timing = [0,12]
base_dir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\20140918\'

figure(3)
clf
hold on
channel = 'RFP'
legend_vec_RFP =  {'WT 4uM 1-NM-PP1',
    'AS 4uM 1-NM-PP1', 
    'AS_2 4uM 1-NM-PP1',
    'AS 0.2 uM 1-NM-PP1' , 
    'AS_2 0.2 uM 1-NM-PP1', 
    'WT SDC',
    'AS SDC', 
    'AS_2 SDC'}
    
fname_save = '20140918_processed_data_SC_2Color.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [0,1,0;  %Green
    0,0,1;  %Blue
    0,0,1;  %Blue
    0,0.5,1;  %Light Blue
    0,0.5,1;  %Light Blue
    0,0,0;  %Black
    0,0,0;  %Black
    0,0,0;  %Black
    ];

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [7,2,4,1];
N_RFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,70,1,8])



figure(4)
clf
hold on
channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP;

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [7,2,4,1];
N_YFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_YFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, KL.MSN2-YFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,70,1,4])



%% GD + Osmo Figure 5, 6, 7, 8
%Data from 20141031_GD_sorb
phases =  {'Pre','Post'} %,'Post'} 
shift_timing = [0,7.5]   
base_dir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\20141031\'
%all locations had 4 sites

%A4 SC 42 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol
%B4 SC 42 t1: SDC, t2: 0.5M Sorbitol
%C4 SC 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol
%D4 SC 42 t1: 4uM 1-NM-PP1, t2: SDC
%E4 SC 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol
%F4 KL 39 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol
%G4 KL 39 t1: SDC, t2: 0.5M Sorbitol 

%second perturbation added at 22.5min after addition of first perturbation
%which was 9.5 min, so 32min absolute time. 

legend_vec = {'42 t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol',
    '42 t1: SDC, t2: 0.5M Sorbitol', 
    '42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
    '42 t1: 4uM 1-NM-PP1, t2: SDC', 
    '37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}

legend_vec = {'KL t1: Gluc DO t2: Gluc DO + 0.5M Sorbitol','KL t1: SDC, t2: 0.5M Sorbitol'}
wellvecSP.SC = {'A4','B4','C4','D4', 'E4'};  %RFP_only 'C8','D8','E8'
wellvecSP.KL = {'F4','G4'};
Nsites = 4;


figure(5)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'t1: Gluc Drop t2: Gluc Drop + 0.5M Sorbitol',
    '42 t1: SDC, t2: 0.5M Sorbitol', 
    '42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
    '42 t1: 4uM 1-NM-PP1, t2: SDC', 
    '37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
    
fname_save = '20141031_processed_data_SC_2Color.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [0.5,0.2,1;  %Purple
    1,0,0;  %Red
    0,1,0;  %Green
    0.5,1,0;  %light Green
    0.5,0.5,0.5;  %Gray
    ];

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2];
N_RFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,90,1,8])



figure(6)
clf
hold on
channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP;

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2];
N_YFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_YFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, KL.MSN2-YFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,90,1,8])


%K.Lactis Plots


figure(7)
clf
hold on
channel = 'RFP'
legend_vec_RFP =  {'t1: Gluc Drop t2: Gluc Drop + 0.5M Sorbitol',
    't1: SDC, t2: 0.5M Sorbitol'}

    
fname_save = '20141031_processed_data_KL_2Color.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [0.5,0.2,1;  %Purple
    1,0,0;  %Red
    ];

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2];
N_RFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('K. Lac cells, KL.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,90,1.2,3.5])



figure(8)
clf
hold on
channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP;

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2];
N_YFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_YFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('K. Lac cells, SC.MSN2-YFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,70,1.2,3.5])



%% AS + Osmo  Figures 9 and 10

%Data from 20141009_AS_sorb
phases =  {'Pre','Post','Post_p2'} %,'Post'} 
shift_timing = [0,9.5,69]  
base_dir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\20141009\'
%all locations had 4 sites

%A3 42 4uM 1-NM-PP1 at time 0
%B3 42 4uM 1-NM-PP1 at time 0 followed by drug + 0.5M Sorbitol @ time 20
%C3 42 4uM 1-NM-PP1 at time 0 followed by drug + 0.25M Sorbitol
%D3 37 4uM 1-NM-PP1 at time 0 followed by drug + 0.5M Sorbitol
%E3 37 4uM 1-NM-PP1 at time 0 followed by drug + 0.25M Sorbitol
%F3 42 0.5M sorbitol at time 0
%G3 42 0.5M sorbitol at time 0 followed by 4uM 1-NM-PP1 @ time 20 

%second perturbation added at 22.5min after addition of first perturbation
%which was 9.5 min, so 32min absolute time. 

legend_vec = {'42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1','42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', '42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol', '37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', '37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol', '42 t1: 0.5M Sorbitol t2: 0.5M Sorbitol', '42 t1: 0.5M Sorbitol t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
wellvecSP.SC = {'A3','B3','C3','D3', 'E3', 'F3', 'G3'};  %RFP_only 'C8','D8','E8'
wellvecSP.KL = {'F8'};
Nsites = 4;


figure(9)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'AS t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1',
    'AS t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol',
    'AS t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol', 
    'WT t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
    'WT t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.25M Sorbitol', 
    'AS t1: 0.5M Sorbitol t2: 0.5M Sorbitol', 
    'AS t1: 0.5M Sorbitol t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}

    
fname_save = '20141009_processed_data_SC_2Color.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [0,0,1;  %Blue
 0.5,0.2,1;  %Purple
 0.8,0.2,1;  %Lighter Purple
 0,.5,0;  %Green
 0,1,0;  %Light Green
 1,0,0;  %Red
 1,0.5,0;  %Orange
    ];

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,3,2,5,4];
N_RFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,90,1,8])



figure(10)
clf
hold on
channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP;

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,3,2,5,4];
N_YFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_YFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, KL.MSN2-YFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,90,1,8])


%% Figure 11 - 14 osmo then glucose dropout


%Data from 20150429_osmo_then_GD_2
phases =  {'Post'} %,'Post'} 
shift_timing = [0]  
base_dir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\20150429_osmo_then_GD_2\'
%all locations had 4 sites
%Gluc Dropout always 0.34% Glucose final, osmo balanced with Sorbitol.
%SC 42: 
%A7 SDC t12 GD
%B7 0.25M sorb t12 GD
%C7 SDC  t24 GD
%D7 0.25M sorb t24 GD
%E7 SDC  t48 GD
%F7 0.25M sorb  t48 GD
%G7 0.125M sorb t12 GD
%H7 0.0625M sorb t12 GD

wellvecSP.SC = {'A7','B7','C7','D7','E7','F7','G7','H7'}; 
wellvecSP.KL = {};
Nsites = 4;

figure(11)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'t0: SDC t36 Gluc DO(0.34%)', 
    't0: 0.25M sorb t36 Gluc DO(0.34%)', 
    't0: SDC t48 Gluc DO(0.34%)', 
    't0: 0.25M sorb t48 Gluc DO(0.34%)', 
    't0: SDC t60 Gluc DO(0.34%)', 
    't0: 0.25M sorb t60 Gluc DO(0.34%)', 
    't0: 0.125M sorb t36 Gluc DO(0.34%)', 
    't0: 0.0625M sorb t36 Gluc DO(0.34%)'
    }
    
fname_save = '20150429_processed_data_SC.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [0,0,1;  %Blue
 0.5,0.2,1;  %Purple
0,0,1;  %Blue
 0.5,0.2,1;  %Purple
 0,0,1;  %Blue
 0.5,0.2,1;  %Purple
 0.65,0.2,1; %int Purple
 0.8,0.2,1;  %Lighter Purple
 ];

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,7,8];
N_RFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, SC.MSN2-RFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,100,1,8])



figure(12)
clf
hold on
channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP;

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,7,8];
N_YFP = length(perm);
legend_vec = legend_vec_RFP(perm);
cmap = cmap_YFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec),1);
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

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('S. Cer cells, KL.MSN2-YFP')
xlabel('time')
ylabel('Nuclear Localization')
axis([0,100,1,4])



return


%% Single Cell Plots
%normalize by transforming all points linearly from the range spanning the
%minimum to the maximum value of the mean for all experiments to 0-1.
norm_val_max.(channels_to_image{1}) = 1;
norm_val_max.(channels_to_image{2}) = 1;
norm_val_min.(channels_to_image{1}) = 100;
norm_val_min.(channels_to_image{2}) = 100;

for si = 1:length(sites);
    site = sites{si};
    site_ind = find(strcmp(site_list,site));
    all_tracks = all_tracks_vec{site_ind};
    all_times = all_times_vec{site_ind};

    for ph = 1:length(phases)
        figure(ph)
        clf
        phase = phases{ph};
        Ntracks = length(all_tracks.(phase));
        Ntimes = length(all_times.(phase));
        times_ind = 1:Ntimes;
        for ch = 1:length(channels_to_image)
            channel = channels_to_image{ch};
            %Build data storage matrix
            sing_cell_tracks_nfmat = zeros(Ntracks,Ntimes);
             %Go through each track and place it in the correct row at the
             %appropriate timepoint.
            clear tracks
            for jj = 1:Ntracks;
                  tracks_row_nf = [all_tracks.(phase)(jj).nf.(channel)];
                  tracks_row_nmi = [all_tracks.(phase)(jj).nmi.(channel)];
                  tracks_row_times = [all_tracks.(phase)(jj).times];
                  sing_cell_tracks_nfmat(jj,[all_tracks.(phase)(jj).times])= tracks_row_nf;
                  tracks(jj).nf = tracks_row_nf;
                  tracks(jj).nmi = tracks_row_nmi;
                  tracks(jj).times = tracks_row_times;
            end
            
            [mean_nf, std_nf] = nf_calcs(times_ind,tracks,[]);
            sing_cell_tracks.(site).(phase).nf_mat.(channel) = sing_cell_tracks_nfmat;
            sing_cell_tracks.(site).(phase).nf_mean.(channel) = mean_nf;
            norm_val_max.(channel) = max(norm_val_max.(channel),max(mean_nf));
            norm_val_min.(channel) = min(norm_val_min.(channel),min(mean_nf));
            sing_cell_tracks.(site).(phase).nf_std.(channel) = std_nf;
            [mean_nmi, std_nmi] = nmi_calcs(times_ind,[]);
            sing_cell_tracks.(site).(phase).nmi_mean.(channel) = mean_nmi;
            sing_cell_tracks.(site).(phase).nmi_std.(channel) = std_nmi;            
            sing_cell_tracks.(site).(phase).nmi_std.(channel) = std_nmi;
            sing_cell_tracks.(site).(phase).times = all_times.(phase);
        end
        
    end
    
    
end

%Normalize by norm_val
%plot
phase = 'Post'
site = 'A8'
%condition = '2% Glu -> 0.5M Sorbitol';
condition = '2% Glu -> no gluc, 0.11M Sorb';

figure(1)
channel = channels_to_image{1};
sing_cell_mat_1 = sing_cell_tracks.(site).(phase).nf_mat.(channel);
sing_cell_mat_1_norm = (sing_cell_mat_1-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));

%sort rows: by max value in the first channel, first 4 timepoints.
[sing_cell_mat_1_norm_sorted,sort_ind] = sortrows(sing_cell_mat_1_norm,[-1,-2,-3,-4]);
zero_shift = (0-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
sing_cell_mat_1_norm_sorted(sing_cell_mat_1_norm_sorted==zero_shift) = nan;

[nr,nc] = size(sing_cell_mat_1_norm_sorted);
pcolor([sing_cell_mat_1_norm_sorted nan(nr,1); nan(1,nc+1)]);
colormap(cool)
shading flat;
set(gca, 'ydir', 'reverse');
colorbar
title(['SC: SC.MSN2 nuclear localization. ', condition, '.'])
ylabel('Cells')
xlabel('Time')
set(gca,'XTick',1:length([sing_cell_tracks.(site).(phase).times])+0.5)
set(gca,'XTickLabel',sprintf('%0.0f|',[sing_cell_tracks.(site).(phase).times]))
caxis([-0.4,4])

figure(2)
channel = channels_to_image{2};
sing_cell_mat_2 = sing_cell_tracks.(site).(phase).nf_mat.(channel);
sing_cell_mat_2_norm = (sing_cell_mat_2-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
%sort using index from previous channel
sing_cell_mat_2_norm_sorted = sing_cell_mat_2_norm(sort_ind,:);
zero_shift = (0-norm_val_min.(channel))/(norm_val_max.(channel)-norm_val_min.(channel));
sing_cell_mat_2_norm_sorted(sing_cell_mat_2_norm_sorted==zero_shift) = nan;

[nr,nc] = size(sing_cell_mat_2_norm_sorted);
pcolor([sing_cell_mat_2_norm_sorted nan(nr,1); nan(1,nc+1)]);
colormap(cool)
shading flat;
set(gca, 'ydir', 'reverse');
colorbar
title(['SC: KL.MSN2 nuclear localization. ', condition, '.'])
ylabel('Cells')
xlabel('Time')
set(gca,'XTick',1:length([sing_cell_tracks.(site).(phase).times])+0.5)
set(gca,'XTickLabel',sprintf('%0.0f|',[sing_cell_tracks.(site).(phase).times]))
caxis([-0.4,4])

figure(3)
clf 
hold on
cmap = cool(length([sing_cell_tracks.(site).(phase).times]));
[Ncells,Ntimes] = size(sing_cell_mat_1);
for jj = 1:Ntimes;
    cell_vec_1 = sing_cell_mat_1(:,jj);
    cell_vec_1 = cell_vec_1(sing_cell_mat_1(:,jj)>0);
    cell_vec_2 = sing_cell_mat_2(:,jj);
    cell_vec_2 = cell_vec_2(sing_cell_mat_1(:,jj)>0); %used same filter intentionally - although should be the same
    scatter(cell_vec_1,cell_vec_2,5,cmap(jj,:),'filled')
    corr_vec(jj) = corr(cell_vec_1,cell_vec_2);
end
xlabel('SC.MSN2(RFP)')
ylabel('KL.MSN2(YFP)')
title(['SC: ',condition])
axis([1,14,1,9])
%scatterplot of data altering color

figure(4)
plot([sing_cell_tracks.(site).(phase).times],corr_vec,'LineWidth',3)
title(['SC: ',condition,'. Correlation SC.MSN2 n.loc. to KL.MSN2 n.loc'])
xlabel('Time')
ylabel('Corr')
axis([10,70,0,1])


%% Old plots

%set colormap (i.e. map = cool(8)) perhaps make use of ColorOrder
%cmap = jet(length(legend_vec));
cmap = [1,0,0;
0,0,0];



figure(1)
clf 
hold on
for jj = 1:length(legend_vec)
    all_tracks = all_tracks_vec{jj};
    all_times = all_times_vec{jj};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5);
        set(p,'Parent',plt_grp(jj))
    end
end



hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Condition')
title('MSN2 Nuclear Localization after media change')


figure(2)
clf 
hold on
for jj = 1:length(legend_vec)
    all_tracks = all_tracks_vec{jj};
    all_times = all_times_vec{jj};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0,'nmi','linewidth',1.5);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Condition')
title('MSN2 median intensity')


figure(3)
clf 
hold on
for jj = 1:length(legend_vec)
    all_tracks = all_tracks_vec{jj};
    all_times = all_times_vec{jj};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_ncells(times,tracks,color_val);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec); %,'Location','NE');
htitle = get(hleg,'Title');
set(htitle,'String','Condition')
title('Number of cells identified')



figure(4)
clf 
hold on
legend_vec_RFP = {'(KL.MSN2) 39y 1 No gluc, 0.111M sorb';
    '(KL.MSN2) 39y 1 0.5M sorb';
    %'(KL.MSN2) 39y 2 0.5M sorb';
    %'(KL.MSN2) 39y 2 No gluc, 0.111M sorb';
    '(KL.MSN2) 29y 1 No gluc, 0.111M sorb';
    '(KL.MSN2) 29y 1 0.5M sorb'}

fname_save = '20140703_processed_data_KL_RFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')   
%all_times_vec_RFP = all_times_vec;
%all_tracks_vec_RFP = all_tracks_vec;
%posvec_RFP = posvec;

% cmap = [0,0,1; 
% 0,1,0.5; 
% 0,.3333,.8333;
% 0,0.6667,0.6667;
% 1,0,0;
% 0,0,0]

cmap = [0,0,1; 
 .5,0.25,0; 
 %0,.3333,.8333;
 %0,0.6667,0.6667;
 0,0,0;
 1,0,0];

perm = [2,1,6,5];

for jj = 1:length(legend_vec_RFP)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5);
        set(p,'Parent',plt_grp(jj))
    end
end


% legend_vec_YFP = {'(SC.MSN2) 39y 1 0.5M sorb';
%     '(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
%     '(SC.MSN2) 39y 2 0.5M sorb';
%     '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'};

legend_vec_YFP = {'(SC.MSN2) 39y 1 No gluc, 0.111M sorb';
    '(SC.MSN2) 39y 1 0.5M sorb';
%    '(SC.MSN2) 39y 2 0.5M sorb';
%    '(SC.MSN2) 39y 2 No gluc, 0.111M sorb'
     };

cmap = [0,0,1; 
 .5,0.25,0; 
 %0,.3333,.8333;
 %0,0.6667,0.6667;
 %0,0,0;
 %1,0,0
 ];

perm = [2,1];
 
fname_save = '20140703_processed_data_KL_YFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
%all_times_vec_YFP = all_times_vec;
%all_tracks_vec_YFP = all_tracks_vec;
%posvec_YFP = posvec;

for jj = 1:length(legend_vec_YFP)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(jj,:);
    kk = jj+length(legend_vec_RFP); %step up plot group
    plt_grp(kk) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        times = all_times.(phases{ph});
        p = plot_meanvalues(times,tracks,color_val,0,'nf','linewidth',1.5,'LineStyle',':');
        set(p,'Parent',plt_grp(kk))
    end
end

%combine legend vectors
legend_vec = [legend_vec_RFP;legend_vec_YFP];
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.MSN2 and SC.MSN2 response to glucose dropout and sorbitol in KL')
xlabel('time')
ylabel('Nuclear Localization')



%For a particular condition
site = 1;   %this is the index in the saved data for A11,KL 39y1 SDC to SDC+0.5M Sorb 
%For a particular phase
phase = 'Post'

%load data for each color
fname_save = '20140703_processed_data_KL_RFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
all_times_vec_RFP = all_times_vec;
all_tracks_vec_RFP = all_tracks_vec;

all_times_RFP =  all_times_vec_RFP{site}.(phase);
all_tracks_RFP = all_tracks_vec_RFP{site}.(phase);

fname_save = '20140703_processed_data_KL_YFP.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec') 
all_times_vec_YFP = all_times_vec;
all_tracks_vec_YFP = all_tracks_vec;

all_times_YFP = all_times_vec_YFP{site}.(phase);
all_tracks_YFP = all_tracks_vec_YFP{site}.(phase);

if length(all_times_YFP) ~= length(all_times_RFP)
    'Error - different number of time points in each channel'
end

times_ind = 1:length(all_times_vec_YFP);

%Find mean 
[nf_mean_RFP, nf_std_RFP] = nf_calcs(times_ind,tracks);
[nf_mean_RFP, nf_std_RFP] = nf_calcs(times_ind,tracks);




profile off
end



