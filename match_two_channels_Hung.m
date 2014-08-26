function celldata = match_two_channels(celldata_1,channel_1,celldata_2,channel_2,L)
%finds a correspondance between points for images of two different
%channels


%Note: celldata locations are sorted in ascending x coordinate order so no
%need to do that here
points_1 = [[celldata_1.Cxloc]',[celldata_1.Cyloc]'];
points_2 = [[celldata_2.Cxloc]',[celldata_2.Cyloc]'];

%Find candidate points 

JJ = size(points_1,1);
KK = size(points_2,1);
dist_mat = zeros(JJ,KK);

for jj = 1:JJ
    for kk = 1:KK
        dist_mat(jj,kk) = sqrt(sum((points_1(jj,:)-points_2(kk,:)).^2));
    end
end

dist_mat_mask = dist_mat < L;
dist_mat(~dist_mat_mask) = Inf;

[match_mat,cost] = Hungarian(dist_mat);

[cells1,cells2] = find(match_mat);
perm_out = [cells1,cells2] ;
perm_out = sortrows(perm_out,1) ;

%for celldata_2_out, use perm_out to reorder correctly (don't try to index
%zero rows)
celldata_2_out = celldata_2(perm_out((perm_out(:,2)>0),2));
%remove zero rows from celldata_1_out
celldata_1_out = celldata_1(perm_out((perm_out(:,1)>0),1));

comb_Cxloc = ([celldata_1_out.Cxloc] + [celldata_2_out.Cxloc])/2;
comb_Cyloc = ([celldata_1_out.Cyloc] + [celldata_2_out.Cyloc])/2;

celldata = struct();
for jj = 1:length(celldata_1_out);
    celldata_fields = fields(celldata_1_out);
    for kk = 1:length(celldata_fields)
        celldata(jj).(celldata_fields{kk}).(channel_1) = celldata_1_out(jj).(celldata_fields{kk});
        celldata(jj).(celldata_fields{kk}).(channel_2) = celldata_2_out(jj).(celldata_fields{kk});
    end
    celldata(jj).Cxloc.Combined = comb_Cxloc(jj);
    celldata(jj).Cyloc.Combined = comb_Cyloc(jj);
end

end
