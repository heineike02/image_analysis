function [fig_out, timevals, mean_val, std_val] = plot_meanvalues(timevals,tracks,channel,color_val,std_flag,val_to_plot, varargin)


if isrow(timevals)
    timevals = timevals';
    'timevals should be column vector'
end
    
%make input parser to parse varargin arguments
pp = inputParser
addOptional(pp,'norm_val',1,@isnumeric)
addOptional(pp,'plot_params',[])
parse(pp,varargin{:});

norm_val = pp.Results.norm_val;
plot_params = pp.Results.plot_params;

times_ind = 1:length(timevals);
if strcmp(val_to_plot,'nf')
    [mean_val, std_val] = nf_calcs(times_ind,tracks,channel);
elseif strcmp(val_to_plot,'nmi')
    [mean_val, std_val] = nmi_calcs(times_ind,tracks,channel);
end
% should std_val be normalized by some other value?
% mean_val and std_val should be column vectors. 

if ~iscolumn(mean_val)
   error('Error. \n mean_val should be a column vector')
end
   
if ~iscolumn(mean_val)
   error('Error. \n std_val should be a column vector')
end
    

mean_val = mean_val./norm_val;

%if std_flag = 1, plot std, 
%if std_flag = 2, plot error bars
%otherwise, just plot mean
if std_flag == 1
    alpha = 0.2;
    fill([timevals ; flipud(timevals)],[mean_val+std_val ; flipud(mean_val-std_val)],color_val, 'FaceAlpha', alpha,'linestyle','none');
    hold on
    fig_out = plot(timevals,mean_val,'Color',color_val,plot_params{:}); % change color or linewidth to adjust mean line
elseif std_flag == 2
    fig_out = errorbar(timevals,mean_val,std_val,'Color',color_val,plot_params{:}); 
else
    fig_out = plot(timevals,mean_val,'Color',color_val,plot_params{:}); % change color or linewidth to adjust mean line
end



%xlabel('time(min)')
%ylabel('Mean Nuclear Localization')

end
    
