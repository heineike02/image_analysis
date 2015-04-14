fname_save = '20150123_processed_data_SC_Msn2.mat';
base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150123_Hog_MSN2_GD_Osmo\';
load([base_dir,fname_save])
pos_list = [all_tracks_vec{3}.Post.pos];
p1 = find(pos_list == 2);  %because of python numbering v.s. matlab numbering p1 from data is labeled as 2 in matlab
tracks_Post_p1 = all_tracks_vec{3}.Post(p1)

track_vec = [16,26,28,33,39,45]

Cxloc_vec = {tracks_Post_p1.Cxloc}
Cyloc_vec = {tracks_Post_p1.Cyloc}
times_vec = {tracks_Post_p1.times}
nf_vec = {tracks_Post_p1.nf}

RFP_img_t1 = imread('C:\Users\Ben\Documents\Data\PKA_project\20150123_Hog_MSN2_GD_osmo\Post\F9-Site_1\img_000000000_RFP_001.tif');

figure(1)
clf 
imagesc(RFP_img_t1)
colormap(gray)
hold on
cmap = jet(length(track_vec));
for jj = 1:length(track_vec)
    xx = Cxloc_vec{track_vec(jj)};
    yy = Cyloc_vec{track_vec(jj)};
    plot(yy,xx,'Color',cmap(jj,:),'LineWidth',2)
end

figure(2)
clf 
hold on
cmap = jet(length(track_vec));
for jj = 1:length(track_vec)
    tt = times_vec{track_vec(jj)}
    nf = nf_vec{track_vec(jj)}.RFP
    plot(tt,nf,'Color',cmap(jj,:),'LineWidth',3)
end
xlabel('Frame')
ylabel('RFP Nuclear Localization')



%% Just Gluc Drop with error bars - not normalized
figure(5)
clf 
hold on

phases =  {'Pre','Post'}

channel = 'RFP'
%legend_vec_RFP = {'(KL.MSN2) t1: Gluc DO t2: Gluc DO + 0.5M Sorb';'(KL.MSN2) t1: SDC, t2: 0.5M Sorb'}
legend_vec_Msn2_RFP = {'SC MSN2-RFP 2% Gluc + .5M Gluc',
    'SC MSN2-RFP 2% Gluc + .5M Sorb',
    'SC MSN2-RFP .11M Sorb',
    'SC MSN2-RFP SDC'}
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol', 
%    '(SC.MSN2) 42 t1: 4uM 1-NM-PP1, t2: SDC', 
%    '(SC.MSN2) 37 t1: 4uM 1-NM-PP1, t2: 4uM 1-NM-PP1 + 0.5M Sorbitol'}
fname_save = '20150123_processed_data_SC_Msn2.mat';
load([base_dir,fname_save],'all_timevals_vec','all_tracks_vec','posvec')
%all_times_vec_Msn2 = all_times_vec;
%all_tracks_vec_Msn2 = all_tracks_vec;
%posvec_Msn2 = posvec;
cmap = [1,0,0;  %Red
 0,1,0;  %Green
 0,0,1;  %Blue
 0,0,0   %Black
 ];  

perm = [3];

%index of normalizing condition
norm_ind = 4
plot_params = {'linewidth',1.5,'LineStyle','-'}

for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap(perm(jj),:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,1,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end


%add Legend vector
legend_vec = legend_vec_Msn2_RFP(perm);
hleg = legend(plt_grp,legend_vec) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SC.MSN2 response in S.Cerevisiae cells')
xlabel('time')
ylabel('Nuclear Localization')
%}



