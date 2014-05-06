function plot_individual_values(tracks)

hold all
for nn = 1:length(tracks)
   nf_vec = [tracks(nn).nf];
   t_vec = [tracks(nn).times];
   plot(t_vec,nf_vec)
end
xlabel('time(min)')
ylabel('Nuclear Localization')
