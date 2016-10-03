Nsites = 6
wellvec = {'A1','B1','C1','D1','E1','F1','A2','B2','C2','D2','E2','F2','G2'} %,'B2','C2','D2','E2','F2','G2'}
%wellvec = {'A5','B5','C5','D5','E5','F5','G5','H5','A6'}

Nwells = length(wellvec);
for jj = 1:Nwells;
    for kk = 1:Nsites;
        posvec{jj,kk} = [wellvec{jj},'-Site_',num2str(kk-1)];
    end
end