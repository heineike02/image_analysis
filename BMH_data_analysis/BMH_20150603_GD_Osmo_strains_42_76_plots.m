function all_tracks = BMH_20150603_GD_Osmo_strains_42_76_plots()
%

profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

base_dir = 'C:\Users\Ben\Documents\Data\PKA_project\20150603_GD_Osmo_strains_42_76\'
%input information from experiment here
species = 'SC' %if I wanted to cycle would make a cell {'SC'}  %Right now not properly cycling through each species

phases =  {'Pre','Post'}%,'Post_p2'} %,'Post'} 
shift_timing = [0,12]    

%std thresh for calling a peak to find a cell
%0.16 seemed to work best for k.lactis on 16JUL images
std_threshSP.KL = 0.16;
std_threshSP.SC = 0.2;

%Tracking parameters
%estimated maximum number of tracks
maxdisp_1x = 4;

% optical amplification 
% 1x or 1.5x
op_amp =  '1.5x'
storeim = 1;


% Msn2 (Strain 42) wells has RFP and YFP
%wellvecSP.SC = {'A6','B6','E6','F6'}; 
%
%Removed A6 site1 
%posvecSP.SC{1,2} = 'NA'


% Hog1 (strain 76) wells: Only has RFP channel
%wellvecSP.SC = {'C6','D6','G6','H6'}; 
%

%Removed C6 site3 
%posvecSP.SC{1,4} = 'NA'

%Removed G6 site2 
%posvecSP.SC{3,3} = 'NA'




%all locations had 4 sites
%Initial Glucose dropout 0.08. 
%osmo balanced with sorbitol 
% Strain 42 has SC.MSN2-YFP and KL.MSN2-RFP
% Strain 76 has Hog1-RFP

%A6	Strain 42 t0  GD 0.25M Sorb t21
%B5	Strain 42 t0 SDC 0.25M Sorb t21
%C5	Strain 76 t0  GD 0.25M Sorb t21
%D5	Strain 76 t0 SDC 0.25M Sorb t21
%E5	Strain 42 t0  GD 0.25M Sorb t36
%F5	Strain 42 t0 SDC 0.25M Sorb t36
%G5	Strain 76 t0  GD 0.25M Sorb t36
%H5	Strain 76 t0 SDC 0.25M Sorb t36


%{
%% Mean time traces

%SC.Msn2-RFP
figure(1)
clf
hold on
channel = 'RFP'
legend_vec_RFP = {'t0: 0.8% Gluc t21: 0.25M Sorb', 
    't0: 2% Gluc t21: 0.25M Sorb', 
    't0: 0.8% Gluc t36: 0.25M Sorb', 
    't0: 2% Gluc t361: 0.25M Sorb'}

fname_save = '20150603_processed_data_SC_MSN2.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [ 0.5,0.2,1;
    0,0,1;  %Blue
    0.8,0.2,1;  %Lighter Purple
    0,0.5,1
 ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3,4] ;
N_RFP = length(perm)
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_RFP),1);
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

hleg = legend(plt_grp,legend_vec_RFP) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SC.Msn2-RFP, Osmo Shock (0.25M Sorbitol) following glucose dropout (2% -> 0.08% osmo balanced w/Sorbitol)')
xlabel('time')
ylabel('Nuclear Localization')

%KL.MSN2-YFP
figure(2)
clf 
hold on
channel = 'YFP'
legend_vec_YFP = legend_vec_RFP
cmap_YFP = cmap_RFP
%perm same as above
legend_vec_YFP = legend_vec_YFP(perm);
cmap = cmap_YFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle',':'};
plt_grp = zeros(length(legend_vec_YFP),1);
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


hleg = legend(plt_grp,legend_vec_YFP) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('KL.Msn2-RFP, Osmo Shock (0.25M Sorbitol) following glucose dropout (2% -> 0.08% osmo balanced w/Sorbitol)')
xlabel('time')
ylabel('Nuclear Localization')

%}

%SC.Hog1-RFP
figure(3)
clf
hold on
channel = []
legend_vec_RFP = {'t0: 0.8% Gluc t21: 0.25M Sorb', 
    't0: 2% Gluc t21: 0.25M Sorb', 
    't0: 0.8% Gluc t36: 0.25M Sorb', 
    't0: 2% Gluc t361: 0.25M Sorb'}

fname_save = '20150603_processed_data_SC_HOG1.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap_RFP = [ 0.5,0.2,1;
    0,0,1;  %Blue
    0.8,0.2,1;  %Lighter Purple
    0,0.5,1
 ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3,4] ;
N_RFP = length(perm)
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_RFP),1);
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

hleg = legend(plt_grp,legend_vec_RFP) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SC.Hog1-RFP, Osmo Shock (0.25M Sorbitol) following glucose dropout (2% -> 0.08% osmo balanced w/Sorbitol)')
xlabel('time')
ylabel('Nuclear Localization')

%plot individual cells for hog
figure(4)
clf
hold on
channel = []
title_vec = {'t0: 0.8% Gluc t21: 0.25M Sorb', 
    't0: 2% Gluc t21: 0.25M Sorb', 
    't0: 0.8% Gluc t36: 0.25M Sorb', 
    't0: 2% Gluc t361: 0.25M Sorb'}

%already loaded
%fname_save = '20150603_processed_data_SC_HOG1.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

perm = [1,2,3,4] ;
N_RFP = length(perm)
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);

for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    subplot(2,2,jj)
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_individual_values(timevals,tracks,channel,'nf');
    end
    title(title_vec{jj})
    axis([0,125,1,6])
end

%filter out hog bright cells

%other ideas for a filter: 
% -  Use the positions I gathered by hand and filter those out.  
% -  Take 2 Standard deviations above normal NF and filter cells that are
% above that most/all of the time
%
% Could try experiment again. 

figure(5)
clf 
% Check Position D6 site 0 v.s. D6 site 1
all_tracks = all_tracks_vec{2};
track_pos_vec = [all_tracks.Post.pos];

subplot(2,2,1)
tracks_D6_site_0 = all_tracks.Post(track_pos_vec == 1);
timevals = all_times.Post;
p = plot_individual_values(timevals,tracks_D6_site_0,channel,'nf');
subplot(2,2,2)
p = plot_individual_values(timevals,tracks_D6_site_0,channel,'nmi');

subplot(2,2,3)
tracks_D6_site_1 = all_tracks.Post(track_pos_vec == 2);
timevals = all_times.Post;
p = plot_individual_values(timevals,tracks_D6_site_1,channel,'nf');
subplot(2,2,4)
p = plot_individual_values(timevals,tracks_D6_site_1,channel,'nmi');

%any plot that has NF above 3 looks like it might be an MSN2 cell

%Verify this by looking at locations of all filtered points: 
nf_thresh = 3;
filter_index = [];
for jj = 1:length(tracks_D6_site_1)
    curr_track_nf = tracks_D6_site_1(jj).nf;
    if sum(curr_track_nf > nf_thresh) > 0
        filter_index = [filter_index, jj];
    end
end
filter_index

%got consitutively bright cell and MSN2 only. 
% Cycle through all locations, positions, and filter out bad tracks (defined as NF above threshold)
% Note that this filter only applies to this dataset as described above.. 

filter_data_file = fopen([base_dir, 'filtered_data.txt'],'w');

nf_thresh = 3;

well_list = {'C6','D6','G6','H6'}

all_tracks_vec_filtered = {}

for jj = 1:length(all_tracks_vec)
    all_tracks = all_tracks_vec{jj};
    all_tracks_filtered = [];
    track_positions = [all_tracks.Post.pos];
    processed_positions = unique(track_positions);
    for kk = 1:length(processed_positions)
        site_tracks = all_tracks.Post(track_positions == processed_positions(kk));
        filter_index = [];
        for nn = 1:length(site_tracks)
            curr_track_nf = site_tracks(nn).nf;
            if sum(curr_track_nf > nf_thresh) > 0  %if any measurements in any track are above 3, add it to the filter index
                filter_index = [filter_index, nn];
            end
        end
        
        fprintf(filter_data_file,'%12s \n',[well_list{jj},'-Site-',int2str(processed_positions(kk))-1]);
        for nn = 1:length(filter_index)
            xloc = mean(site_tracks(filter_index(nn)).Cxloc);
            yloc = mean(site_tracks(filter_index(nn)).Cyloc);
            fprintf(filter_data_file,'xloc %6.2f yloc %6.2f \n',xloc, yloc);
        end
        %remove wells
        site_tracks_filtered = site_tracks;
        site_tracks_filtered(filter_index) = [];
        all_tracks_filtered = [all_tracks_filtered, site_tracks_filtered];
    end
    all_tracks_vec{jj}.Post_filtered = all_tracks_filtered;
    all_times_vec{jj}.Post_filtered = all_times_vec{jj}.Post;
end

fclose(filter_data_file);

%SC.Hog1-RFP filtered
figure(5)
phases = {'Pre','Post_filtered'};
clf 
hold on
channel = []
legend_vec_RFP = {'t0: 0.8% Gluc t21: 0.25M Sorb', 
    't0: 2% Gluc t21: 0.25M Sorb', 
    't0: 0.8% Gluc t36: 0.25M Sorb', 
    't0: 2% Gluc t361: 0.25M Sorb'}
cmap_RFP = [ 0.5,0.2,1;
    0,0,1;  %Blue
    0.8,0.2,1;  %Lighter Purple
    0,0.5,1
 ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2,3,4] ;
N_RFP = length(perm)
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_RFP),1);
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

hleg = legend(plt_grp,legend_vec_RFP) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('SC.Hog1-RFP, Osmo Shock (0.25M Sorbitol) following glucose dropout (2% -> 0.08% osmo balanced w/Sorbitol) \n filtered')
xlabel('time')
ylabel('Nuclear Localization')

%plot individual cells for hog
figure(6)
clf
hold on
channel = []
title_vec = {'t0: 0.8% Gluc t21: 0.25M Sorb', 
    't0: 2% Gluc t21: 0.25M Sorb', 
    't0: 0.8% Gluc t36: 0.25M Sorb', 
    't0: 2% Gluc t361: 0.25M Sorb'}

%already loaded
%fname_save = '20150603_processed_data_SC_HOG1.mat';
%load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

perm = [1,2,3,4] ;
N_RFP = length(perm)
legend_vec_RFP = legend_vec_RFP(perm);
cmap = cmap_RFP(perm,:);

for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    subplot(2,2,jj)
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_individual_values(timevals,tracks,channel,'nf');
    end
    title(title_vec{jj})
    axis([0,125,1,3])
end
suptitle('Filtered')


%Display 


%{
%C6 Site 0: 
1	351	33
2	371	36
3	130	352
4	140	364
5	322	34

C6-Site-0 
xloc 351.35 yloc 115.60 
xloc 358.89 yloc 130.94


%C6 Site 1
1	365	47
2	382	46
3	313	263 
4	334	264
5	283	271
6	292	231

   C6-Site-1 
xloc  63.45 yloc 353.32 
xloc 259.79 yloc 318.89 
xloc  54.83 yloc 349.83 
xloc  43.13 yloc 352.18 

%C6 Site 2
1	395	241
2	410	230
3	421	221
4	401	259

   C6-Site-2 
xloc 230.40 yloc 412.40 
xloc 216.23 yloc 423.04 

%C6 Site 3 - blob, removed from analysis

%D6 Site 0 - no YFP Cells

%D6 Site 1
1	148	501
2	162	493 - 5
3	244	455
4	252	473 - 1,3

2 x = 400, y = 86
4 x = 400, y = 86
6 x = 394, y = 86
This was a bright red splotch.  

xloc 474.67 yloc 263.83 
xloc 400.29 yloc  85.29 
xloc 472.67 yloc 254.06 
xloc 396.32 yloc  75.74 
xloc 496.69 yloc 160.45 
xloc 393.94 yloc  84.78 

Criteria seems sufficient to remove wells. 

%D6 Site 2 - no YFP Cells

%D6 Site 3
1	447	187
2	464	190
3	485	190
4	492	174

xloc 194.00 yloc 455.13 
xloc 174.11 yloc 498.56 
xloc 185.25 yloc 488.75 
xloc 188.55 yloc 487.58 
xloc 187.83 yloc 450.79 

%G6 Site 0 
1	332	10
2	348	13
3	338	17

   G6-Site-0 
xloc  18.22 yloc 334.11 
xloc  15.71 yloc 332.71 

%G6 Site 1 - no YFP cells

%G6 Site 3
1	326	293
2	327	311
3	312	303
4	315	279

Site 3 removed


%G6 Site 3 - last two move before sticking
1	273	123
2	291	115
3	263	170
4	280	181

   G6-Site-3 
xloc 120.00 yloc 298.13 

None

%H6 Site 0 no YFP cells

%H6 Site 1  half of one YFP cell on edge

   H6-Site-1 
xloc 111.69 yloc 474.69 

%H6 Site 2
1	477	149
2	460	160
3	37	503
4	49	498
5	18	502

   H6-Site-2 
xloc 156.15 yloc 481.25 
xloc 159.33 yloc 476.00 
xloc 157.38 yloc 487.31 
xloc 154.37 yloc 487.63 
xloc 135.15 yloc 484.69 
xloc 158.72 yloc 464.00

%H6 Site 3
1	414	187
2	428	206
3	437	185
4	452	170
5	417	220
%bright thing
6	398	30
7	391	43
8	407	42
9	389	52
10	382	43

   H6-Site-3 
xloc 189.03 yloc 443.84 
xloc 207.83 yloc 422.33 
xloc 213.60 yloc 421.20 
xloc 190.27 yloc 418.60 
xloc 198.67 yloc 427.00 
xloc 204.95 yloc 429.38 
xloc  59.64 yloc 406.21 
xloc  55.75 yloc 405.63 

%}




