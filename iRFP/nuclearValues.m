function [nuclearVals] = nuclearValues(outMask,imgRFP, showFig)

% apply mask to grayscale image
imgRFP1 = outMask.*double(imgRFP);
% obtain intensity values for image
ff = bwconncomp(imgRFP1);
stats2= regionprops(ff,imgRFP1, 'Area','Centroid','PixelList', 'PixelValues');
% sum of pixels (norm to area)
sumFF=[];
% median of pixels (norm to area)
medFF=[];
% mean of pixels (norm to area)
meanFF=[];
for i = 1:length([stats2.Area])
    temp1 = stats2(i).PixelValues;
    temp2 = stats2(i).Area;
    sumFF(i) = sum(temp1);
    medFF(i) = median(temp1);
    meanFF(i) = mean(temp1);
end
% show figs
if showFig
    figure(11); imagesc(imgRFP1);
end

nuclearVals.sum = sumFF;
nuclearVals.med = medFF;
nuclearVals.mean = meanFF;

end