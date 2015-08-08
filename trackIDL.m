function tracks = trackIDL(timecoursedata, track_params)

time_inds = track_params.time_inds;
maxdisp = track_params.maxdisp;
track_memory = track_params.track_memory;
min_points_for_traj = track_params.min_points_for_traj;
maxcells = track_params.maxcells;
ipdir = track_params.ipdir;

%adds IDL Track function to path
%Windows
path([ipdir,'IDL_Particle_Tracking'],path);

%See track documentation for more info on these parameters
param.mem = track_memory;

param.dim = 2;
param.good = min_points_for_traj;

%allows error messages to print
param.quiet = 0;




nFrames = length(timecoursedata);

%calculates number of channels from number of fields in the nf field of the
%first cell in the first frame of timecoursedata. 
channels = fields(timecoursedata(1).celldata(1).nf);
nChan = length(channels);


%initialize position vector
xyzs = zeros(maxcells*nFrames, 3+2*nChan);

mm = 1; 
for jj = 1:length(timecoursedata)
    celldata = timecoursedata(jj).celldata;
    nCells = length(celldata);
    for kk = 1:nCells
        xyzs(mm,1) = celldata(kk).Cxloc.Combined;
        xyzs(mm,2) = celldata(kk).Cyloc.Combined;
        for nn = 1:nChan
            xyzs(mm,1+2*nn) = celldata(kk).nf.(channels{nn});
            xyzs(mm,2+2*nn) = celldata(kk).nmi.(channels{nn});
        end
        xyzs(mm,3+2*nChan) = time_inds(jj);
        mm = mm + 1;
    end
end


%Get rid of additional zeros at the end of the vector
xyzs = xyzs(sum(xyzs,2)>0,:);

% %plot xy positions
% figure(1)
% clf
% hold on
% colorvec = jet(length(timecoursedata));
% for jj = 1:length(timecoursedata)
%     xy = xyzs(xyzs(:,5)==jj,[1,2]);
%     plot(xy(:,1),xy(:,2),'Marker','x','LineStyle','none','Color',colorvec(jj,:))
% end

track_raw = track(xyzs,maxdisp,param);
track_inds = track_raw(:,end);

nTracks = track_inds(end);

for jj = 1:nTracks
    track_inds_jj = (track_inds ==jj);
    tracks(jj).Cxloc = track_raw(track_inds_jj,1);
    tracks(jj).Cyloc = track_raw(track_inds_jj,2);
    for kk = 1:nChan
        tracks(jj).nf.(channels{kk}) = track_raw(track_inds_jj,1+2*kk);
        tracks(jj).nmi.(channels{kk}) = track_raw(track_inds_jj,2+2*kk);
    end
    tracks(jj).times = track_raw(track_inds_jj,3+2*nChan);
    tracks(jj).length = sum(track_inds_jj);
end

% xyzs = zeros(maxcells*nFrames, 5);
% 
% mm = 1; 
% for jj = 1:length(timecoursedata)
%     celldata = timecoursedata(jj).celldata;
%     nCells = length(celldata);
%     for kk = 1:nCells
%         xyzs(mm,1) = celldata(kk).Cxloc;
%         xyzs(mm,2) = celldata(kk).Cyloc;
%         xyzs(mm,3) = celldata(kk).nf;
%         xyzs(mm,4) = celldata(kk).nmi;
%         xyzs(mm,5) = time_inds(jj);
%         mm = mm + 1;
%     end
% end
% 
% %Get rid of additional zeros at the end of the vector
% xyzs = xyzs(sum(xyzs,2)>0,:);
% 
% %plot xy positions
% % figure(1)
% % clf
% % hold on
% % colorvec = jet(length(times));
% % for jj = 1:length(times)
% %     xy = xyzs(xyzs(:,5)==jj,[1,2]);
% %     plot(xy(:,1),xy(:,2),'Marker','x','LineStyle','none','Color',colorvec(jj,:))
% % end
% 
% track_raw = track(xyzs,maxdisp,param);
% track_inds = track_raw(:,end);
% 
% nTracks = track_inds(end);
% 
% for jj = 1:nTracks
%     track_inds_jj = (track_inds ==jj);
%     tracks(jj).Cxloc = track_raw(track_inds_jj,1);
%     tracks(jj).Cyloc = track_raw(track_inds_jj,2);
%     tracks(jj).nf = track_raw(track_inds_jj,3);
%     tracks(jj).nmi = track_raw(track_inds_jj,4);
%     tracks(jj).times = track_raw(track_inds_jj,5);
%     tracks(jj).length = sum(track_inds_jj);
% end
% 





end
