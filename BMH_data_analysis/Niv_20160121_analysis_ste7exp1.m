function [] = BMH_20151215_KL_NMPP1_Osmo()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)




%%
base_dir = 'C:\Users\Ben\Documents\Data\Niv\22_01_16Ste7exp1\'

%phases =  {'Post'}%,'Post_p2'} %,'Post'} 
%shift_timing = [0]    



%Well	Strain             Condition            Pert1                   Pert2 plate
%A4     78(KL TPK2/3 AS)   NMPP1 + 0.5M Sorb    100ul SDC + 8uM NMPP1	200ul 1M Sorb + 4uM NMPP1
%B4     78(KL TPK2/3 AS)   0.5M Sorb only       100ul SDC               200ul 1M Sorb
%C4     78(KL TPK2/3 AS)   NMPP1 + 0.25M Sorb	100ul SDC + 8uM NMPP1	200ul 0.5M Sorb + 4uM NMPP1
%D4     78(KL TPK2/3 AS)   0.25M Sorb           100ul SDC               200ul 0.5M Sorb
%E4     74(KL WT)          NMPP1 + 0.5M Sorb	100ul SDC + 8uM NMPP1	200ul 1M Sorb + 4uM NMPP1
%F4     74(KL WT)          0.5M Sorb Only       100ul SDC               200ul 1M Sorb
%G4     52(SC TPK1/2/3 AS) NMPP1 + 0.5M Sorb	100ul SDC + 8uM NMPP1	200ul 1M Sorb + 4uM NMPP1
%H4     52(SC TPK1/2/3 AS) 0.5M Sorb only       100ul SDC               200ul 1M Sorb


fname_save = 'processed_data_CD_1to3.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

titles = {'C','D'}
t = 1
for jj = [1:2]
    figure(1)
    subplot(1,2,jj)
    Ncells = length(all_tracks_vec{jj}.Exp)
    nf_vec = []
    nmi_vec = []
    for kk = 1:Ncells
        times = all_tracks_vec{jj}.Exp(kk).times
        if ismember(t,times)
            t_ind = find(times == t);
            nf_vec = [nf_vec, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_ind)];
            nmi_vec = [nmi_vec, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_ind)];
        end
    end
    plot(nf_vec,nmi_vec,'x')
    xlabel('nf vec')
    ylabel('nmi vec')
    title(titles{jj})
    axis([0 16 0 11])
end

all_tracks_vec;
all_times_vec;



fname_save = 'processed_data_CD_allpos.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

titles = {'C','D'}
t = 1
for jj = [1:2]
    Ncells = length(all_tracks_vec{jj}.Exp)
    nf_vec = [];
    nmi_vec = [];
    for kk = 1:Ncells
        times = all_tracks_vec{jj}.Exp(kk).times
        if ismember(t,times)
            t_ind = find(times == t);
            nf_vec = [nf_vec, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_ind)];
            nmi_vec = [nmi_vec, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_ind)];
        end
    end
    figure(2)
    subplot(1,2,jj)
    [X,bins] = hist3([nmi_vec',nf_vec'],[20,20]);
    imagesc(X)
    xlabel('nf vec')
    ylabel('nmi vec')
    title(titles{jj})
        
    
    figure(3)
    subplot(1,2,jj)
    plot(nf_vec,nmi_vec,'x')
    xlabel('nf vec')
    ylabel('nmi vec')
    title(titles{jj})
    axis([0 20 0 20])
    
    
    figure(4)
    subplot(1,2,jj)
    hist(nf_vec, 30)
    title('NF')
    xlim([0,25])
    
    figure(5)
    subplot(1,2,jj)
    hist(nmi_vec, 30)
    title('NMI')
    xlim([0,25])
end





return
%% RFP channel 

channel = 'RFP'
legend_vec_RFP = {'TPK2/3 AS, 4uM 1NMPP1 + 0.5M Sorb','TPK2/3 AS, 0.5M Sorb' ,'TPK2/3 AS, 4uM 1NMPP1 + 0.25M Sorb','TPK2/3 AS, 0.25M Sorb', 'WT, 4uM 1NMPP1 + 0.5M Sorb','WT, 0.5M Sorb'  }
cmap_RFP = [       1,0,1;  %Purple
     1,0,0;  %Red
     0.5,0,1;  %Light Purple
     0.5,0,0;  %Light Red
     0,0,0;  %Black
     0,0,0;  %Black
          ];  

linestyle_vec_RFP = {'-','-','-','-','--','-'}
marker_vec_RFP = {'None','None','None','None','None','None'}

%KL.MSN2 for all strains
figure(1)
clf
hold on

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [2,4,6,1,3,5];
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
axis([0,130,1,3])
title('KL.MSN2-RFP in KL Cells')
xlabel('time')
ylabel('Nuclear Localization')



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



