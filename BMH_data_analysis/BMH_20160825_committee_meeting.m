function [] = BMH_20160825_committee_meeting()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)


%Goal is to plot all the different conditions on one plot

%%Part one: process data for TM experiment

%all_tracks = BMH_20140409_KL_SC_H2O2_data_process_20160825()

base_dir = 'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\20160826_Ctte_Mtg_summary_fig\'

species = {'SCer', 'KLac'}
%listed as condition, filename, index for data

%First Convert all datasets that are in the old style to the new style
%where channel is always included. 

%{
old_style_data = {'processed_data_20140605_gluc_drop_SC','processed_data_20140605_gluc_drop_KL', 'processed_data_20140812_RFP_gal_gly_SC', 'processed_data_20140812_RFP_gal_gly_KL' }
new_fname = {'processed_data_20140605_gluc_drop_SC_new', 'processed_data_20140605_gluc_drop_KL_new', 'processed_data_20140812_RFP_gal_gly_SC_new', 'processed_data_20140812_RFP_gal_gly_KL_new'}
for jj = 1:length(old_style_data)
    fname = [base_dir, old_style_data{jj}];
    load(fname,'all_times_vec','all_tracks_vec','posvec');
    all_tracks_vec = add_channel_label(all_tracks_vec,'RFP');
    save([base_dir,new_fname{jj}],'all_times_vec','all_tracks_vec','posvec')    
end

%}



conditions_SC = {{'0.5 M Sorbitol', 'processed_data_20140605_gluc_drop_SC_new',1},
    {'Gluc ->3% Gly','processed_data_20140812_RFP_gal_gly_SC_new',2},
    {'1.0mM H2O2','processed_data_20140409_H202_SC',3},
    {'Glucose removal','processed_data_20140605_gluc_drop_SC_new',3},
    {'Gluc -> Gal','processed_data_20140812_RFP_gal_gly_SC_new',1},
    {'Gluc -> Gly','processed_data_20140812_RFP_gal_gly_SC_new',3},
    {'Gluc -> Sorb','processed_data_20140605_gluc_drop_SC_new',4},
    {'1NMPP1 (AS strain)','processed_data_20161008_NMPP1_SC',2}};

conditions_KL = {{'0.5 M Sorbitol', 'processed_data_20140605_gluc_drop_KL_new',1},
    {'Gluc ->3% Glycerol','processed_data_20140812_RFP_gal_gly_KL_new',2},
    {'1.0mM H2O2','processed_data_20140409_H202_KL',3},
    {'Gluc removal','processed_data_20140605_gluc_drop_KL_new',3},
    {'Gluc -> Gal','processed_data_20140812_RFP_gal_gly_KL_new',1},
    {'Gluc -> Gly','processed_data_20140812_RFP_gal_gly_KL_new',3},
    {'Gluc -> Sorb','processed_data_20140605_gluc_drop_KL_new',4},
    {'1NMPP1 (AS strain)','processed_data_20151113_NMPP1_KL',2}};

conditions_species = {conditions_SC,conditions_KL};


figure(1)
clf
hold on

%for each species
for sp = 1:length(species)
    conditions = conditions_species{sp};
    %for each condition
    for jj = 1:length(conditions)
        jj
        %open up the file
        fname = [base_dir, conditions{jj}{2}];
        load(fname,'all_times_vec','all_tracks_vec','posvec');
        
        %plot the data
        phases = fieldnames(all_tracks_vec{1});
        all_tracks = all_tracks_vec{conditions{jj}{3}};
        all_times = all_times_vec{conditions{jj}{3}};
        
        %plot_params = {'linewidth',1.5,'LineStyle',linestyle_vec_RFP{perm(jj)},'Marker',marker_vec_RFP{perm(jj)}};
        subplot(2,length(conditions_SC), jj + 8*(sp-1))
        hold on
        jj+8*(sp-1)
        plot_params = {'linewidth', 1.5};
        for ph = 1:length(phases)
            tracks = all_tracks.(phases{ph});
            timevals = all_times.(phases{ph});
            [fig_out, timevals, mean_val, std_val] = plot_meanvalues(timevals,tracks,'RFP','k',0,'nf','plot_params',plot_params);
            %set(p,'Parent',plt_grp(jj))
            if sp==1
                title(conditions{jj}{1})
                axis([0,60,1.5,8])
            elseif sp == 2
                axis([0,60,1.5,3.3])
                xlabel('time(min)')
            end
            if jj == 1
                ylabel('Nuclear Localization')
            end
            
        end
    end 
end

return

% K.Lactis response
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20151215_KL_NMPP1_Osmo\'



phases =  {'Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    



%Well	Strain             Condition            Pert1                   Pert2 plate
%A4     78(KL TPK2/3 AS)   NMPP1 + 0.5M Sorb    100ul SDC + 8uM NMPP1	200ul 1M Sorb + 4uM NMPP1
%B4     78(KL TPK2/3 AS)   0.5M Sorb only       100ul SDC               200ul 1M Sorb
%C4     78(KL TPK2/3 AS)   NMPP1 + 0.25M Sorb	100ul SDC + 8uM NMPP1	200ul 0.5M Sorb + 4uM NMPP1
%D4     78(KL TPK2/3 AS)   0.25M Sorb           100ul SDC               200ul 0.5M Sorb
%E4     74(KL WT)          NMPP1 + 0.5M Sorb	100ul SDC + 8uM NMPP1	200ul 1M Sorb + 4uM NMPP1
%F4     74(KL WT)          0.5M Sorb Only       100ul SDC               200ul 1M Sorb
%G4     52(SC TPK1/2/3 AS) NMPP1 + 0.5M Sorb	100ul SDC + 8uM NMPP1	200ul 1M Sorb + 4uM NMPP1
%H4     52(SC TPK1/2/3 AS) 0.5M Sorb only       100ul SDC               200ul 1M Sorb


fname_save = 'processed_data_KL.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

all_tracks_vec;
all_times_vec;

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



