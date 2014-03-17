function tracks = trackIDLu(datainput,times, maxdisp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION - tracks cell
%
%INPUTS:
%           datainput - structure with cell segmentation information
%           times - a vector of time points
%           maxdisp - the maximum allowed gap between two traces in time
%
%OUTPUTS:
%           tracks - a structure of the tracks found
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%See track documentation for more info on these parameters

%number of frames that the particle can disappear for
param.mem = 2; 

%xy data
param.dim = 2; 

%track length
if length(times)<5
    param.good = 1;
else
    param.good = floor(2*length(times)./3); 
end

param.quiet = 0;

%maximum number of cells expected per frame
maxcells = 1000;

%number of frames
nFrames = length(datainput);

%initialize position vector
xyzs = zeros(maxcells*nFrames, 5);

%create the structure for the tracking algorithm
mm = 1; 
%iterate through frames 
for jj = 1:length(datainput) 
    stats = datainput(jj).stats;
    nCells = length(stats); 
    %iterate through each cell
    for kk = 1:nCells 
        %locX
        xyzs(mm,1) = stats(kk).Centroid(1);
        %locY
        xyzs(mm,2) = stats(kk).Centroid(2);
        %just place holder
        xyzs(mm,3) = 1;
        %mean intensity
        xyzs(mm,4) = stats(kk).MeanIntensity;
        %timepoint
        xyzs(mm,5) = times(jj);
        mm = mm + 1;
    end
end
 
%Get rid of additional zeros at the end of the vector
xyzs = xyzs(sum(xyzs,2)>0,:);

%tracking happens here
track_raw = celltrack(xyzs,maxdisp,param);
track_inds = track_raw(:,end);

nTracks = track_inds(end);

%storing tracked information
for jj = 1:nTracks
    track_inds_jj = (track_inds ==jj);
    tracks(jj).locX = track_raw(track_inds_jj,1);
    tracks(jj).locY = track_raw(track_inds_jj,2);
    tracks(jj).nf = track_raw(track_inds_jj,3);
    tracks(jj).nmi = track_raw(track_inds_jj,4);
    tracks(jj).times = track_raw(track_inds_jj,5);
    tracks(jj).length = sum(track_inds_jj);
end