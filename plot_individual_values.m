function fig_out = plot_individual_values(timevals,tracks,channel,val_to_plot,varargin)

hold all
xlabel('time(min)')
for nn = 1:length(tracks)
   if isempty(channel)
      y_vec = [tracks(nn).(val_to_plot)];
   else
      y_vec = [tracks(nn).(val_to_plot).(channel)];
   end
   time_inds = tracks(nn).times;
   t_vec = timevals(time_inds);
   fig_out = plot(t_vec,y_vec);
   %patchline
end


