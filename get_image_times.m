function [time_inds, time_vals] = get_image_times(fname_conv,imdir,channel_to_image,pos_fnames_nn,images,time_calc);
%time_calc:  Tells the program how to calculate each time value.  If the
%input is a number, then that is interpreted as a dt between images.  If
%the input is a filename, then that is interpreted as metadata that gives
%exact time values.  

nTimes = length(images.(channel_to_image));
time_inds = 1:nTimes;
if isnumeric(time_calc)
   %method of calculating time is dt 
   time_vals = time_calc*(0:(nTimes-1));
elseif ischar(time_calc) % either JSO or Micromanager
   time_vals = zeros(size(time_inds));
   if strcmp(fname_conv,'JSO')
       metadata = regexp(importdata(time_calc),'\t','split');
       time_fnames = cell(size(metadata));
       for jj = 1:length(metadata)
            time_fnames{jj} = metadata{jj}{9};
       end
   
       for jj = 1:nTimes % time that corresponds to the t# for 1 well
            jj;
            fname = images(1).(channel_to_image){jj};
            %acqdata file numbers filenames in the metadata differently (start
            %time at 0) than it names the actual filenames (starting at 1)
            %fname_times = char(strcat(channel_to_image,'_',pos_fnames_nn,'_t',int2str(channel_mod-1+(mm-1)*Nchan),'.tiff'));      
            [fname_parts,fname_time] = regexpi(fname,'_t([0-9]*).','split','tokens');
            fname_metadata = [fname_parts{1},'_t',num2str(str2num(fname_time{1}{1})-1),'.',fname_parts{2}];
            metadata_ind = strcmp(time_fnames,fname_metadata);
            time = str2double(metadata{metadata_ind}{1})/60.0; %in minutes
            time_vals(jj) = time;
       end
   elseif strcmp(fname_conv,'Micromanager')
       [~,Channel,~,Frame,Time] =  import_metadata_parsed([imdir,pos_fnames_nn,'\',time_calc]);
       chan_ind = strcmp(Channel,channel_to_image);
       time_vals = Time(chan_ind);
       
       %%Old routine that uses Matlab's JSON parser - very slow!!!
%        %This routine assumes all metadata is in order - if it is out of
%        %order, we would have to check the filename
%        metadata = loadjson([imdir,pos_fnames_nn,'\',time_calc]);
%        %get starting time
%        tin = metadata.Summary.Time;
%        t0 = get_time(tin);
%            
%        field_names = fieldnames(metadata);
%        %Keep only fields that belong to images
%        field_names = field_names(~strcmp(field_names,'Summary'));
%        mm = 1;
%        for jj = 1:length(field_names)
%            if strcmp(metadata.(field_names{jj}).Channel,channel_to_image)
%                tin = metadata.(field_names{jj}).Time;
%                t1 = get_time(tin);
%                time = etime(t1,t0)/60; %gets elapsed time in minutes
%                time_vals(mm) = time;
%                mm = mm+1;
%            end
%        end
    end

end

end

% function tout = get_time(tin)
%     %converts time in micromanager metadata format to matlab format - now
%     %done in python
%     t = zeros(1,6);
%     tout = regexp(tin,' ','split');
%     tout = tout{2};
%     tout = regexp(tout,':','split');
%     for jj = 1:3
%        t(jj+3) = str2double(tout{jj});
%     end
%     tout = t;
% end
