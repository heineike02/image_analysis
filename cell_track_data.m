function [tracks,times,images] = cell_track_data(imdir,channel_to_image,pos_fnames_nn,Nchan,shift_spacing_ph,dt,circ,imbg,siz,storeim,rad,std_thresh,maxdisp,t);

    images = struct();
    files=dir(char(strcat(imdir,channel_to_image,'_',pos_fnames_nn,'_t*.tiff'))); 
    nTimes = size(files,1);
    %Get initial number
    init_fname = files(1).name;
    init_ind = strfind(init_fname,'_t');
    dot_ind = strfind(init_fname,'.');
    init_num = str2num(init_fname(init_ind+2:(dot_ind-1)));
    %picks out the first number and checks what it is modulo the number of channels 
    channel_mod = mod(init_num,Nchan)+Nchan*(mod(init_num,Nchan)==0) ;
    %This is a pain because the numbering is not 01,02,03 and
    %messes up the ordering of the files. 
    for mm = 1:nTimes 
       images(mm).(channel_to_image) = char(strcat(channel_to_image,'_',pos_fnames_nn,'_t',int2str(channel_mod+(mm-1)*Nchan),'.tiff'));
    end
    
    %Shift times for each phase
    t = t + shift_spacing_ph;
    times = dt*[0:nTimes-1]+t;
    t = times(end);
    
    timecoursedata = CovMapCellsBMH(imdir ,images, {channel_to_image}, circ, imbg , siz , storeim, rad, std_thresh, times);
    tracks = trackIDL_BMH(timecoursedata,times, maxdisp);

    %filter const bright cells
    tracks = filterTracksBMH(tracks);

end

