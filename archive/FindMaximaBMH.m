function [maxx,maxy]=FindMaximaBMH(im,pointx,pointy,dis)

maxx = [];
maxy = [];
%sets margin of 10 around border of image, filters out points within edge
%boundaries
edge=10;
a=(find(pointx>edge & pointx<size(im,1)-edge & pointy>edge & pointy<size(im,2)-edge));
pointy=pointy(a);
pointx=pointx(a);

k=0;
for i=1:length(pointx)
    %makes a neighborhood of distance dis around a point
    lipx=im(pointx(i)-dis:pointx(i)+dis,pointy(i)-dis:pointy(i)+dis);
    %sets point to max if the maximum in the neighborhood is at that point
    max_lipx = max(lipx(:));
    if max_lipx-im(pointx(i),pointy(i))==0
        k=k+1;
        maxx(k)=pointx(i);
        maxy(k)=pointy(i);
    end
end

end