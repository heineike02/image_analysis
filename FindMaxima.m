function [maxx,maxy]=FindMaxima(im,pointx,pointy,dis,rm_dupes,edge)

%sets margin around border of image, filters out points within edge
%boundaries
a=(find(pointx>edge & pointx<size(im,1)-edge & pointy>edge & pointy<size(im,2)-edge));
pointy=pointy(a);
pointx=pointx(a);

NN = length(pointx);
maxx = zeros(NN,1);
maxy = zeros(NN,1);

for i=1:NN
    %makes a neighborhood of distance dis around a point
    nbhd = im(pointx(i)-dis:pointx(i)+dis,pointy(i)-dis:pointy(i)+dis);
    %sets maxx and maxy to the max around that point
    [xmax_vec,xmax_ind_vec] = max(nbhd);
    [max_val,ymax_ind] = max(xmax_vec);
    xmax_ind = xmax_ind_vec(ymax_ind);
    
    shiftx = xmax_ind-(dis+1);
    shifty = ymax_ind-(dis+1);
    
    maxx(i)=pointx(i)+shiftx;
    maxy(i)=pointy(i)+shifty;
        
end

if rm_dupes == 1;
    %removes duplicates
    maxxy = [maxx,maxy];
    maxxy = unique(maxxy,'rows');
    maxx = maxxy(:,1);
    maxy = maxxy(:,2);
end

end