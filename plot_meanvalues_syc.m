function [fig_out, timevals, mean_val, std_val] = plot_meanvalues_syc(timevals,tracks,channel,color_val,std_flag,val_to_plot, axisVals, varargin)

%make input parser to parse varargin arguments
pp = inputParser
addOptional(pp,'norm_val',1,@isnumeric)
addOptional(pp,'plot_params',[])
parse(pp,varargin{:});

norm_val = pp.Results.norm_val; 
plot_params = pp.Results.plot_params;

times_ind = 1:length(timevals);
if strcmp(val_to_plot,'nf')
    [mean_val, std_val] = nf_calcs_syc(times_ind,tracks,channel);
elseif strcmp(val_to_plot,'nmi')
    [mean_val, std_val] = nmi_calcs(times_ind,tracks,channel);
end
% should std_val be normalized by some other value?

mean_val = mean_val./norm_val; 

%if std_flag = 1, plot std, otherwise, just plot mean
if std_flag == 1
    alpha = 0.2;
    fill([timevals fliplr(timevals)],[mean_val'+std_val' fliplr(mean_val'-std_val')],color_val, 'FaceAlpha', alpha,'linestyle','none');
    hold on
end
%fig_out = plot(timevals,mean_val,'Color',color_val,plot_params{:}); % change color or linewidth to adjust mean line
%fig_out = errorbar(timevals,mean_val,std_val,'Color',color_val,plot_params{:});
fig_out = shadedErrorBar(timevals, mean_val, std_val);

axis(axisVals)
%xlabel('time(min)')
%ylabel('Mean Nuclear Localization')

end
    
