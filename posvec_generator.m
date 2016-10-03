Nsites = 25
wellvec = {'A3','B3','C3','D3','E3','F3','G3','H3','G2','H2'} %,'B2','C2','D2','E2','F2','G2'}
%wellvec = {'A5','B5','C5','D5','E5','F5','G5','H5','A6'}

Nwells = length(wellvec);
for jj = 1:Nwells;
    for kk = 1:Nsites;
        posvec{jj,kk} = [wellvec{jj},'-Site_',num2str(kk-1)];
    end
end