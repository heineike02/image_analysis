function [cellV,cellx,celly]=AlignCellsKMEANS(cff,Cxloc,Cyloc,tim,ci)

a=find(tim==min(tim));
[IDX, C] = KMEANS([Cxloc',Cyloc'], length(a),'emptyaction','singleton','START',[Cxloc(a)',Cyloc(a)']);

%test for squareness of input
if mod(sqrt(size(cff,2)),1)~=0
    yinc=1;
    siz=(sqrt(size(cff,2)/2));
else
    yinc=0;
    siz=(sqrt(size(cff,2)));
end

for i=1:max(IDX)
    m(i)=length(find(IDX==i));
end

k=0;
fsc=zeros(size(IDX));

for i=1:length(find(m>250));
     ma=find(m>250);
    a=find(IDX==ma(i));
        for j=1:250;
            ab=find(IDX==ma(i) & tim'==j);
            ml(j)=length(ab);
        end
    nc(i)=max(quantile(ml,0.95),1);
    kvs=kmeans([Cxloc(a)',Cyloc(a)'],nc(i),'emptyaction','singleton');
        for j=1:(nc(i))
            mv=length(find(kvs==j));
            if mv>0.75*(max(tim)-min(tim));
                k=k+1; 
                fsc(a(find(kvs==j)))=k;
            end
        end
end

a=find(fsc>0);

if yinc; cellV=nan(k,2,max(tim)); else; cellV=nan(k,max(tim));end

for i=1:length(a)
    nf=NEnrich(reshape(cff(a(i),1:siz^2),siz,siz),ci);
    if yinc
        nfy=NEnrich(reshape(cff(a(i),1+siz^2:2*siz^2),siz,siz),ci);
        cellV(fsc(a(i)),1,tim(a(i)))=nf;
        cellV(fsc(a(i)),2,tim(a(i)))=nfy;
    else
        cellV(fsc(a(i)),tim(a(i)))=nf; 
    end
    cellx(fsc(a(i)),tim(a(i)))=Cxloc(a(i));
    celly(fsc(a(i)),tim(a(i)))=Cyloc(a(i));
end

end


function maxrl=runlength(dat)
maxrl=0;    
k=0;
smrl=0;
for i=1:length(dat)
        if ~isnan(dat(i))
            maxrl=maxrl+1;
        else
            if maxrl>0
                k=k+1;
                smrl(k)=maxrl;
            end
            maxrl=0;
        end
end
k=k+1;
smrl(k)=maxrl;

maxrl=max(smrl);
end

function maxrl=maxgapsize(dat)
maxrl=0;    
k=0;
smrl=0;
for i=1:length(dat)
        if isnan(dat(i))
            maxrl=maxrl+1;
        else
            if maxrl>0
                k=k+1;
                smrl(k)=maxrl;
            end
            maxrl=0;
        end
end
k=k+1;
smrl(k)=maxrl;

maxrl=max(smrl);
end


function nf=NEnrich(im,ci)
    siz=size(im,1);
    %im=conv2(im,ones(3,3),'same');    
    %c2=ci(1:13,2:14);
    c2=ci;
    a=conv2(im,c2,'same');
    
    [in,xy]=max(a(:));
    yloc=(mod(xy-1,siz(1))+1);
    xloc=floor(xy/siz(1))+1;
    
    for i=1:siz;
        matx(i,:)=1:siz;
        maty(:,i)=1:siz;
    end
    
    dis=sqrt((xloc-matx(:)).^2+(yloc-maty(:)).^2);
    ne=sort(im(dis<6),'descend');
    nf=mean(ne(1:5))./median(ne);
    
end

function m=fastnanmean(x)
    m=mean(x(~isnan(x)));
end
