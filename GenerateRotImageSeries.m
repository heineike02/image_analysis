function rf=GenerateRotImageSeries(cf,siz,sizf)

rf=zeros(size(cf,1)*10,2*sizf^2);

for i=1:size(cf,1)
    for j=0:35;
        rif=rotate_image(j*10, reshape(cf(i,1:siz^2),siz,siz));
        yif=rotate_image(j*10, reshape(cf(i,1+siz^2:2*siz^2),siz,siz));
        rif=rif((size(rif,1)-sizf)/2 +1:(size(rif,1)-sizf)/2+sizf,(size(rif,2)-sizf)/2 +1:(size(rif,2)-sizf)/2+sizf);
        yif=yif((size(yif,1)-sizf)/2 +1:(size(yif,1)-sizf)/2+sizf,(size(yif,2)-sizf)/2 +1:(size(yif,2)-sizf)/2+sizf);
        
         rf((i-1)*10+j+1,:) = [rif(:)',yif(:)']';
    end
end