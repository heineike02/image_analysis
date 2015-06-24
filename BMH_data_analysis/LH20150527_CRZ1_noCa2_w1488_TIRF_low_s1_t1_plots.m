function all_tracks = LH20150527_CRZ1_noCa2_w1488_TIRF_low_s1_t1_plots()
%

profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)

base_dir = 'C:\Users\Ben\Documents\Data\0-CRZ1-test\'
%input information from experiment here

phases =  {'Exp'}%,'Post_p2'} %,'Post'} 
shift_timing = [0]    


%Strain BMH42 - SC.MSN2-RFP,  KL.MSN2-YFP 
%wellvecSP.SC = {'A9','B9','C9','D9','E9','F9','G9','H9'};
%

%all locations had 4 sites
% Glucose Dropout (0.1 % from 2%) followed by Osmo shock (0.35M) at various times.
% osmo balanced with sorbitol relative to osmotic stress condition.
%
%
%A8	t9: GD + sorb   
%B8	t9: GD      
%C8	t9: SDC + sorb
%D8	t9: GD  	t15: sorb
%E8	t9: GD      t21: sorb
%F8	t9: GD      t33: sorb
%G8	t9: GD      t45: sorb
%H8	t9: GD      t69: sorb



%% SC.Msn2-RFP Mean Time Traces
figure(1)
clf
hold on
channel = []
legend_vec = {'s1 and s2',
's3 and s4'}

fname_save = '20150527-CRZ1-noCa2_w1488 TIRF low_processed_data.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
cmap = [ 0,0,1;  %Blue
    1,0,0; %Red
    ];  

%cmap = jet(length(legend_vec_RFP));
%perm = 1:length(legend_vec_RFP);
perm = [1,2]  ;
N = length(perm)
legend_vec_plot = legend_vec(perm);
cmap_plot = cmap(perm,:);
plot_params = {'linewidth',1.5,'LineStyle','-'};
plt_grp = zeros(length(legend_vec_plot),1);
for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    color_val = cmap_plot(jj,:);
    plt_grp(jj) = hggroup;
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_meanvalues(timevals,tracks,channel,color_val,0,'nf','plot_params',plot_params);
        set(p,'Parent',plt_grp(jj))
    end
end

hleg = legend(plt_grp,legend_vec_plot) %,'Location','NE');
htitle = get(hleg,'Title');
%set(htitle,'String','Condition')
title('CRZ1-test')
xlabel('time')
ylabel('Nuclear Localization')


%% SC.MSN2-RFP plot individual cells 
figure(2)
clf
hold on
channel = []
title_vec = {'s1 and s2',
's3 and s4'}

%already loaded

perm = [1,2] ;

for jj = 1:length(perm)
    all_tracks = all_tracks_vec{perm(jj)};
    all_times = all_times_vec{perm(jj)};
    subplot(1,2,jj)
    for ph = 1: length(phases)
        tracks = all_tracks.(phases{ph});
        timevals = all_times.(phases{ph});
        p = plot_individual_values(timevals,tracks,channel,'nf');
    end
    axis([0,30,1,2.6])
    title(title_vec{jj})
end


