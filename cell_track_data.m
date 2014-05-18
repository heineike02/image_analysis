function [tracks,times,images] = cell_track_data(imdir,channel_to_image,pos_fnames_nn,Nchan,circ,imbg,siz,storeim,rad,std_thresh,maxdisp,time_calc)

%time_calc:  Tells the program how to calculate each time value.  If the
%input is a number, then that is interpreted as a dt between images.  If
%the input is a filename, then that is interpreted as metadata that gives
%exact time values.  

if isnumeric(time_calc)
    time_calc_method = 'dt'
elseif ischar(time_calc)
    time_calc_method = 'metadata'
end

    images = struct();
    files=dir(char(strcat(imdir,channel_to_image,'_',pos_fnames_nn,'_t*.tiff'))); 
    nTimes = size(files,1);
    times_ind = 1:nTimes;
    if strcmp(time_calc_method,'dt')
        times = dt*(0:(nTimes-1));
    elseif strcmp(time_calc_method,'metadata')
        metadata = regexp(importdata(time_calc),'\t','split');
        time_fnames = cell(size(metadata));
        for jj = 1:length(metadata)
            time_fnames{jj} = metadata{jj}{9};
        end
    end
    
    %Get initial number
    init_fname = files(1).name;
    init_ind = strfind(init_fname,'_t');
    dot_ind = strfind(init_fname,'.');
    init_num = str2double(init_fname(init_ind+2:(dot_ind-1)));
    %picks out the first number and checks what it is modulo the number of channels 
    channel_mod = mod(init_num,Nchan)+Nchan*(mod(init_num,Nchan)==0) ;
    %This is a pain because the numbering is not 01,02,03 and
    %messes up the ordering of the files. 
    for mm = 1:nTimes 
       fname = char(strcat(channel_to_image,'_',pos_fnames_nn,'_t',int2str(channel_mod+(mm-1)*Nchan),'.tiff'));
       images(mm).(channel_to_image) = fname;
       if strcmp(time_calc_method,'metadata')
           %acqdata file numbers filenames in the metadata differently (start
           %time at 0) than it names the actual filenames (starting at 1)
           fname_times = char(strcat(channel_to_image,'_',pos_fnames_nn,'_t',int2str(channel_mod-1+(mm-1)*Nchan),'.tiff'));      
           metadata_ind = find(strcmp(time_fnames,fname_times));
           time = str2double(metadata{metadata_ind}{1})/60.0; %in minutes
           times(mm) = time;
       end
    end
    
    timecoursedata = CovMapCellsBMH(imdir ,images, {channel_to_image}, circ, imbg , siz , storeim, rad, std_thresh, times_ind);
    tracks = trackIDL_BMH(timecoursedata,times_ind, maxdisp);
    %note: tracks use index of times rather than actual times
    
    %filter const bright cells
    tracks = filterTracksBMH(tracks);

end

