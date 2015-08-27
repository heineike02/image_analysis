Nsites = 4
wellvec = {'A11','B11','C11','D11','E11','F11','G11','H11','C10'} %,'B2','C2','D2','E2','F2','G2'}
%wellvec = {'A5','B5','C5','D5','E5','F5','G5','H5','A6'}

Nwells = length(wellvec);
for jj = 1:Nwells;
    for kk = 1:Nsites;
        posvec{jj,kk} = [wellvec{jj},'-Site_',num2str(kk-1)];
    end
end