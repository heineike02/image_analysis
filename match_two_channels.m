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
%set up matrix of possible matches
matches_per_row = sum(dist_mat_mask,2);
match_mat = zeros(JJ,2,max(matches_per_row));

for jj = 1:JJ
    nposs = matches_per_row(jj);
    if nposs == 0
       ['no matching cell for celldata_1 cell ',int2str(jj)]
    else
       row_matches = find(dist_mat_mask(jj,:)>0);
       for kk = 1:nposs;
          match_mat(jj,:,kk) = [jj,row_matches(kk)]; 
       end
    end
end

%build permutations and find distances
NP = prod(matches_per_row(matches_per_row>1));
permutations = zeros(JJ,2,NP);

if NP == 1
    permutations(:,:,1) = match_mat(:,:,1);
else

    mult_match = find(matches_per_row>1);
    perm_dims = cell(1,length(mult_match));
    for mm = 1:length(mult_match);
        perm_dims{mm} = 1:matches_per_row(mult_match(mm));
    end
    %perm_grid = cell(matches_per_row(mult_match)');
    perm_grid = cell(1,numel(perm_dims));
    [perm_grid{:}] = ndgrid(perm_dims{:});

    %number of elements in perm_grid should be NP
    if length(perm_grid{1}(:)) ~= NP
        'Error: permutation number is incorrect'
    end

    for np = 1:NP
        permutations(:,:,np) = match_mat(:,:,1);
        for mm = 1:length(mult_match);
            permutations(mult_match(mm),:,np) = match_mat(mult_match(mm),:,perm_grid{mm}(np));
        end
    end

end
%for each permutation calculate the distance
perm_dist = inf*ones(NP,1);
for np = 1:NP
    perm = permutations(:,:,np);
    perm_cols = perm(:,2);
    %determine if there is a conflict in the column index
    unique_perm_cols = unique(perm_cols);
    count_perm_cols = hist(perm_cols,unique_perm_cols);
    dupe_ind = find(count_perm_cols > 1 );
    if not(isempty(dupe_ind))   
        %if so keep the pair with the lowest distance, replace the thrown out row with 0,0.
        for kk = 1:length(dupe_ind)
            dupe_kk_ind = find(perm_cols == unique_perm_cols(dupe_ind(kk)));
            for dd = 1:length(dupe_kk_ind);
                if (perm(dupe_kk_ind(dd),1)==0) & (perm(dupe_kk_ind(dd),2)==0)
                    d_options(dd) = L;
                else
                    d_options(dd) = dist_mat(perm(dupe_kk_ind(dd),1),perm(dupe_kk_ind(dd),2));
                end
            end
            [~,min_d_ind]= min(d_options);
            %replace all but minimum distance duplicate with (0,0)
            dupe_kk_ind(min_d_ind) = [];
            perm(dupe_kk_ind,:) = zeros(length(dupe_kk_ind),2);
        end
    end
    
    %Add up the distances for all rows
    %if (0,0) is in any row, assign a distance of length L
    perm_d = 0;
    for jj = 1:JJ
        if (perm(jj,1)==0) && (perm(jj,2)==0)
            perm_d = perm_d + L;
        else
            perm_d = perm_d + dist_mat(perm(jj,1),perm(jj,2));
        end
    end
    perm_dist(np) = perm_d; 
    permutations(:,:,np) = perm;
end

[~,min_perm_ind] = min(perm_dist);
perm_out = permutations(:,:,min_perm_ind);

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
