function [acqd,strain,numw,numc]=readAcqDat(fi,sfi)
%open acqdat.txt
fid=fopen(fi);
i=0;
while ~feof(fid)
   i=i+1;
   l=fgetl(fid);
   
   c=(regexpi(l,'\t'));
   el=length(l);
   
   t(i)=str2num(l(1:c(1)-1));
   
   k=1;
   wellnum(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   LocX(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   LocY(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   LocZ(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   PFoff(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   ch(i)=cellstr(l(c(k)+1:c(k+1)-1)); k=k+1;
   exp(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
   
   if length(c)>9
       ledsetval(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
       ledoutval(i)=str2num(l(c(k)+1:c(k+1)-1)); k=k+1;
       fname(i)=cellstr(l(c(k)+1:el(1))); 
   else
       fname(i)=cellstr(l(c(k)+1:el(1)));  
   end
end

fclose(fid);

numw=length(unique(wellnum));
numc=length(unique(ch));    
acqd.wellnum=wellnum;
acqd.ch=ch;
acqd.time=t;
acqd.exp=exp;
if length(c)>9
    acqd.LEDset=ledsetval;
    acqd.LEDout=ledoutval;
end
acqd.fname=fname;

i=0;
if ~isempty(sfi)
    %open strains.txt
    fid=fopen(sfi);
    while ~feof(fid)
        i=i+1;
        strain(i)=cellstr(fgetl(fid));
    end
    fclose(fid)
else
    strain(1:(numw))=cellstr(numw);
end

end