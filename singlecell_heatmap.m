function fig_out = singlecell_heatmap(channel,normval,all_tracks, all_times)

%if normval is 'auto' the normalization value is automatically calculated. 
if strcmp(normval,'auto')
    %calculate normalization value automatically. 
end

fig_out = figure(1)
for ph = 1:length(phases)
    

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

        [mean_nf, std_nf] = nf_calcs(times_ind,tracks);
        sing_cell_tracks.(site).(phase).nf_mat.(channel) = sing_cell_tracks_nfmat;
        sing_cell_tracks.(site).(phase).nf_mean.(channel) = mean_nf;
        norm_val_max.(channel) = max(norm_val_max.(channel),max(mean_nf));
        norm_val_min.(channel) = min(norm_val_min.(channel),min(mean_nf));
        sing_cell_tracks.(site).(phase).nf_std.(channel) = std_nf;
        [mean_nmi, std_nmi] = nmi_calcs(times_ind,tracks);
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
site = 'B12'
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
title(['KL: KL.MSN2 nuclear localization. ', condition, '.'])
ylabel('Cells')
xlabel('Time')
set(gca,'XTickLabel',sprintf('%0.2f|',[sing_cell_tracks.(site).(phase).times]))
caxis([-0.3,3.1])
