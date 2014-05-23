function tracks = cell_track_data(imdir,images,times_ind,times_vals,channel_to_image,circ,imbg,siz,storeim,rad,std_thresh,maxdisp)
    timecoursedata = CovMapCellsBMH(imdir ,images, channel_to_image, circ, imbg , siz , storeim, rad, std_thresh, times_vals);
    tracks = trackIDL_BMH(timecoursedata,times_ind, maxdisp);
    %note: tracks use index of times rather than actual times
    %filter const bright cells
    tracks = filterTracksBMH(tracks);
end