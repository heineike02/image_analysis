%For a given set of all_tracks and all_times, gives the number of cells for
%each desired time

%Scer
x_shift = 7.5;

%Klac
%x_shift = 8.2;

field_name = 'Exp' %Klac: 'Post'  %Scer: 'Exp';

max_time = 13.0;

%SCer
legend_vec_RFP = {'PKA-AS, No Drug','WT, No Drug', 'AS 10nM', 'AS 25 nM', 'AS 100nM', 'WT 100nM', 'AS 250nM', 'AS 1uM', 'PKA-AS, 1-NM-PP1', 'WT, 1-NM-PP1', 'AS 10uM', 'WT 10uM','WT DMSO'};

%KLac
%legend_vec_RFP = {'TPK2-AS, 1-NM-PP1', 'TPK2/3-AS +Cas9, 1-NM-PP1','PKA-AS, 1-NM-PP1','TPK2-AS SDC','TPK2/3-AS +Cas9, SDC','PKA-AS, No Drug','WT, 1-NM-PP1', 'WT, No Drug'}

%Scer
perm = [9,10,1,2];

%Klac
%perm = [3,7,6,8];

for jj = 1:length(perm)
    perm_ind = perm(jj);
    all_times_well = all_times_vec{1,perm_ind}.(field_name);
    all_times_well_subset = all_times_well(all_times_well< x_shift + max_time);
    
    all_tracks_well = all_tracks_vec{1,perm_ind}.(field_name) ;
    
    cells_in_well_times = [];
    for kk = 1:length(all_times_well_subset)
        %for each time, find how many cells are present
        Ncells = 0;
        for mm = 1:length(all_tracks_well)
            track_times = all_tracks_well(mm).times;
            if not(isempty(intersect(kk,track_times)))
                Ncells = Ncells + 1;
            end
        end
        cells_in_well_times = [cells_in_well_times,Ncells];
    end
    
    cells_in_well{jj} = cells_in_well_times;
    
    ['Number of cells in ',legend_vec_RFP{perm_ind}, ' by time']
    cells_in_well{jj}
end
