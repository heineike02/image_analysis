function [cellV,cellx,celly,cfo]=AlignCells(cff,Cxloc,Cyloc,tim,ci,cfi,nuc)

cellx=nan(max(tim),5000);
celly=nan(max(tim),5000);

mitim=min(tim);
matim=max(tim);

if nargin<7; nuc=1; end
if nargin<6; cfi=0; end
if nargin<5 | isempty(ci);
    S=load('circ.mat')
    circ=S.circ;
    ci=circ(1:13,2:14);
end

%test for squareness of input
if mod(sqrt(size(cff,2)),1)~=0
    yinc=1;
    siz=(sqrt(size(cff,2)/2));
else
    yinc=0;
    siz=(sqrt(size(cff,2)));
end

if size(cff,2)==2 | size(cff,2)==4
   storim=0; 
   if size(cff,2)==2
       yinc=0;
   else
       yinc=1;
   end
else
   storim=1;
end


if yinc
    cellV=nan(max(tim),2,500);
    if cfi; cfo=nan(max(tim),500,2*siz.^2); end
else
    cellV=nan(max(tim),1,500);
    if cfi; cfo=nan(max(tim),500,siz.^2); end
end    
    
    
a=find(tim==min(tim)); mnc=length(a);
for j=1:length(a);
    if storim
        [nf,nmi]=NEnrich(reshape(cff(a(j),1:siz^2),siz,siz),ci);
        if yinc; [nfy,nmiy]=NEnrich(reshape(cff(a(j),1+siz^2:2*(siz^2)),siz,siz),ci); end
    else
        nf=cff(a(j),1);
        nmi=cff(a(j),2);
        if yinc; 
            nfy=cff(a(j),3);
            nmiy=cff(a(j),4); 
        end        
    end
    cellx(min(tim),j)=Cxloc(a(j));
    celly(min(tim),j)=Cyloc(a(j));
    if nuc
        cellV(min(tim),1,j)=nf; 
        if yinc; cellV(min(tim),2,j)=nfy; end
    else
        cellV(min(tim),1,j)=nmi;
        if yinc; cellV(min(tim),2,j)=nmiy; end
    end
    
    if cfi; cfo(min(tim),j,:)=(cff(a(j),:)); end
end

    for j=1:mnc
        cavX(j)=nanmean(cellx(:,j));
        cavY(j)=nanmean(celly(:,j));
    end

for i=min(tim)+1:max(tim)
    a=find(tim==i);
    for j=1:length(a);
    if storim
        [nf,nmi]=NEnrich(reshape(cff(a(j),1:siz^2),siz,siz),ci);
        if yinc; [nfy,nmiy]=NEnrich(reshape(cff(a(j),1+siz^2:2*(siz^2)),siz,siz),ci); end
    else
        nf=cff(a(j),1);
        nmi=cff(a(j),2);
        if yinc; 
            nfy=cff(a(j),3);
            nmiy=cff(a(j),4); 
        end  
    end
        dx=cellx(i-1,:)-Cxloc(a(j));
        dy=celly(i-1,:)-Cyloc(a(j));
        dt=sqrt(dx.^2+dy.^2);
        [mv,mi]=min(dt);
        if mv<12
            mif=mi;
        else
            dx=cavX-Cxloc(a(j));
            dy=cavY-Cyloc(a(j));
            dt=sqrt(dx(:).^2+dy(:).^2);
            [mv,mi]=min(dt);
            if mv<8
                mi=floor(mi./(i-1))+1;
                mif=mi;
            else
                mnc=mnc+1;
                mif=mnc;
            end
        end
        cellx(i,mif)=Cxloc(a(j));
        celly(i,mif)=Cyloc(a(j));
        if nuc
            cellV(i,1,mif)=nf; 
            if yinc; cellV(i,2,mif)=nfy; end
        else
            cellV(i,1,mif)=nmi;
            if yinc; cellV(i,2,mif)=nmiy; end
        end

        if cfi; cfo(i,mif,:)=(cff(a(j),:)); end
        
        
    end
    %nanmean(cellx(i,:)-cellx(i-1,:))
    %nanmean(celly(i,:)-celly(i-1,:))
    for j=1:mnc
        cavX(j)=nanmean(cellx(mitim:i,j));
        cavY(j)=nanmean(celly(mitim:i,j));
    end
end

cellx=cellx(:,1:mnc);
celly=celly(:,1:mnc);
cellV=cellV(:,:,1:mnc);
if cfi; cfo=cfo(:,1:mnc,:); end

for i=1:mnc
    l(i)=sum(~isnan(cellx(:,i)));
    maxrl(i)=runlength(cellx(:,i));
    maxgl(i)=maxgapsize(cellx(:,i));
end

a=find(l>(max(tim)-min(tim))/2 & maxrl>(max(tim)-min(tim))/4);


cellx=cellx(:,a);
celly=celly(:,a);
cellV=cellV(:,:,a);
if cfi;  cfo=cfo(:,a,:); else; cfo=[]; end


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


function [nf,nmi]=NEnrich(im,ci)
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
    nmi=median(ne);
end

function m=fastnanmean(x)
    m=mean(x(~isnan(x)));
end
