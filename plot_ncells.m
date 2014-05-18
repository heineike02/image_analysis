function fig_out = plot_ncells(times,tracks,color_val)
%since number of cells may be above 256, I am using uint16 data type
times_ind = 1:length(times);
ncells = zeros(size(times),'uint16');

for jj = 1:length(tracks)
    track_times = tracks(jj).times;
    for kk = 1:length(times_ind);
        ncells(kk) = ncells(kk)+ uint16(sum(track_times == times_ind(kk)));
    end
end
    
fig_out = plot(times,ncells,'Color',color_val,'linewidth',1.5); % change color or linewidth to adjust mean line
  
xlabel('time(min)')
ylabel('Number of cells')

end
    
