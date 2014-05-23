function tracks = trackIDL_BMH(timecoursedata,times, maxdisp);

%adds IDL Track function to path
%Windows
path(path,'C:\Users\Ben\Google Drive\UCSF\ElSamad_Lab\PKA\Image_analysis\Image_processing_scripts\IDL_Particle_Tracking');

%See track documentation for more info on these parameters
param.mem = 2;
param.dim = 2;
if length(times)<5
    param.good = 1;
else
    param.good = 4;
end
param.quiet = 0;

%maximum number of cells expected per frame
maxcells = 200;


nFrames = length(timecoursedata);

%initialize position vector
xyzs = zeros(maxcells*nFrames, 4);

mm = 1; 
for jj = 1:length(timecoursedata)
    celldata = timecoursedata(jj).celldata;
    nCells = length(celldata);
    for kk = 1:nCells
        xyzs(mm,1) = celldata(kk).Cxloc;
        xyzs(mm,2) = celldata(kk).Cyloc;
        xyzs(mm,3) = celldata(kk).nf;
        xyzs(mm,4) = celldata(kk).nmi;
        xyzs(mm,5) = times(jj);
        mm = mm + 1;
    end
end

%Get rid of additional zeros at the end of the vector
xyzs = xyzs(sum(xyzs,2)>0,:);

%plot xy positions
% figure(1)
% clf
% hold on
% colorvec = jet(length(times));
% for jj = 1:length(times)
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
    tracks(jj).nf = track_raw(track_inds_jj,3);
    tracks(jj).nmi = track_raw(track_inds_jj,4);
    tracks(jj).times = track_raw(track_inds_jj,5);
    tracks(jj).length = sum(track_inds_jj);
end




end
