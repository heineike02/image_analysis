function fig_out = plot_meanvalues(times,tracks,color_val,std_flag,val_to_plot)

times_ind = 1:length(times);
if strcmp(val_to_plot,'nf')
    [mean_val, std_val] = nf_calcs(times_ind,tracks);
elseif strcmp(val_to_plot,'nmi')
    [mean_val, std_val] = nmi_calcs(times_ind,tracks);
end

%if std_flag = 1, plot std, otherwise, just plot mean
if std_flag == 1
    alpha = 0.2;
    fill([times fliplr(times)],[mean_val'+std_val' fliplr(mean_val'-std_val')],color_val, 'FaceAlpha', alpha,'linestyle','none');
    hold on
end
fig_out = plot(times,mean_val,'Color',color_val,'linewidth',1.5); % change color or linewidth to adjust mean line
  
xlabel('time(min)')
ylabel('Mean Nuclear Localization')

end
    
