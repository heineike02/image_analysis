function [maxx,maxy]=FindLocPoint(im,pointx,pointy,dis,sc)

edge=max(10,dis);
a=(find(pointx>edge & pointx<size(im,1)-edge & pointy>edge & pointy<size(im,2)-edge));
pointy=pointy(a);
pointx=pointx(a);

k=0;
for i=1:length(pointx)
    lipx=im(pointx(i)-dis:pointx(i)+dis,pointy(i)-dis:pointy(i)+dis);
    lstd=std(lipx);
    lmean=mean(lipx);
    if im(pointx(i),pointy(i))-lmean-sc*lstd>0
        k=k+1;
        maxx(k)=pointx(i);
        maxy(k)=pointy(i);
    end
end

end