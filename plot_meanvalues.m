function fig_out = plot_meanvalues(times,tracks,color_val,std_flag)

%if std_flag = 1, plot std, otherwise, just plot mean
[mean_nf, std_nf] = nf_calcs(times,tracks);

if std_flag == 1
    alpha = 0.2;
    fill([times fliplr(times)],[mean_nf'+std_nf' fliplr(mean_nf'-std_nf')],color_val, 'FaceAlpha', alpha,'linestyle','none');
    hold on
end
fig_out = plot(times,mean_nf,'Color',color_val,'linewidth',1.5); % change color or linewidth to adjust mean line
  
xlabel('time(min)')
ylabel('Mean Nuclear Localization')

end
    
