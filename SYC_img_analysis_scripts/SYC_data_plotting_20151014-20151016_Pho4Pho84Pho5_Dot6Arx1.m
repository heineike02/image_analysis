% SYC data plotting 20151014-20151016 Quick Version
%%
%% 20151016 RFP
load('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20151014-20151016_Pho4Pho84Pho5_Dot6Arx1_Lights\20151016_Pho4Pho5_Dot6Arx1_NOlight\BlueLightNO_RFP_20151014.mat')
A=all_tracks_vec{1}.BlueLightNo;
B=all_tracks_vec{2}.BlueLightNo;
lengthA = length(A);
lengthB = length(B);
imFreq = 3;
% no blue light
figure(1); subplot(2,1,1); for i=1:lengthA; plot(imFreq.*(A(i).times-1),A(i).nf.RFP); hold on; end; xlabel('time (min)'); ylabel('nuc/cyt'); xlim([0 90]); legend('Pho4 Pho5'); title('no blue light');
subplot(2,1,2); for i=1:lengthB; plot(imFreq.*(B(i).times-1),B(i).nf.RFP); hold on; end; xlabel('time (min)'); ylabel('nuc/cyt'); xlim([0 90]); legend('Dot6 Arx1');

load('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20151014-20151016_Pho4Pho84Pho5_Dot6Arx1_Lights\20151016_Pho4Pho5_Dot6Arx1_NOlight\Dot6Arx1_Pho4Pho5_BL_RFP_20151016.mat')
A=all_tracks_vec{1}.Pho4Pho5;
B=all_tracks_vec{1}.Dot6Arx1;
lengthA = length(A);
lengthB = length(B);
imFreq = 3;
% 45min blue light at 25au
figure(1); subplot(2,1,1); for i=1:lengthA; plot(imFreq.*(A(i).times-1),A(i).nf.RFP); hold on; end; xlabel('time (min)'); ylabel('nuc/cyt'); xlim([0 90]); legend('Pho4 Pho5'); title('45min 25au BL');
subplot(2,1,2); for i=1:lengthB; plot(imFreq.*(B(i).times-1),B(i).nf.RFP); hold on; end; xlabel('time (min)'); ylabel('nuc/cyt'); xlim([0 90]); legend('Dot6 Arx1');

%%
%% 20151015 RFP
load('C:\Users\susanychen\Google Drive\UCSF Graduate Research (El-Samad)\DECODING PKA BY DOWNSTREAM EFFECTORS\20151014-20151016_Pho4Pho84Pho5_Dot6Arx1_Lights\20151015_Pho4Pho84_BL\Pho4Pho84_BL_RFP_20151015.mat')
A=all_tracks_vec{1}.BlueLightNo;
B=all_tracks_vec{1}.BlueLight21;
C=all_tracks_vec{1}.BlueLight45;
lengthA = length(A);
lengthB = length(B);
lengthC = length(C);
imFreq = 3;
% Light at different lengths @ 25au
figure(1); subplot(3,1,1); for i=1:lengthA; plot(imFreq.*(A(i).times-1),A(i).nf.RFP); hold on; end; ylabel('nuc/cyt'); xlim([0 90]); ylim([1 4]); title('Pho4 Pho84: No BL'); %xlabel('time (min)');
subplot(3,1,2); for i=1:lengthB; plot(imFreq.*(B(i).times-1),B(i).nf.RFP); hold on; end; xlabel('time (min)'); ylabel('nuc/cyt'); xlim([0 90]); ylim([1 4]); title('Pho4 Pho84: 21min @ 25au BL');
subplot(3,1,3); for i=1:lengthC; plot(imFreq.*(C(i).times-1),C(i).nf.RFP); hold on; end; xlabel('time (min)'); ylabel('nuc/cyt'); xlim([0 90]); ylim([1 4]); title('Pho4 Pho84: 45min @ 25au BL');

%%
%% 20151016 YFP
cd \\elsamad.ucsf.edu\data\instrumentation\microscope\SYC\20151016_Pho4Pho5_Dot6Arx1_NOlight\Pho4Pho5_BlueLight45min10inten25au\
% extract YFP files
files = dir('*YFP*');
nTimes = size(files,1);
% extract values
A=[];
for i = 1:length(files)
    [mat,tok,ext]=regexpi(files(i).name,'\d+','match', 'tokens');
    A(i) = str2num(mat{2});
end
[valA, indA] = sort(A);
filesO = files(indA); % ordered

ratioIm = [];
for i = 1:length(filesO);
    im = imread(filesO(i).name);
    im = double(im);
    imMed = median(im(:));
    imStd = std(im(:));
    imThresh = imMed+1.*imStd;
    indThresh = find(im>imThresh);
    oneStdIm = median(im(indThresh));
    ratioIm(i) = oneStdIm/imMed;
end

numPts = 30;
figure(2); plot(3.*(([1:numPts]-1)),ratioIm); hold on;

cd \\elsamad.ucsf.edu\data\instrumentation\microscope\SYC\20151016_Pho4Pho5_Dot6Arx1_NOlight\BlueLightNo\
% extract YFP files
files = dir('YFP_p1*');
nTimes = size(files,1);
% extract values
A=[];
for i = 1:length(files)
    [mat,tok,ext]=regexpi(files(i).name,'\d+','match', 'tokens');
    A(i) = str2num(mat{2});
end
[valA, indA] = sort(A);
filesO = files(indA); % ordered

ratioIm = [];
for i = 1:length(filesO);
    im = imread(filesO(i).name);
    im = double(im);
    imMed = median(im(:));
    imStd = std(im(:));
    imThresh = imMed+1.*imStd;
    indThresh = find(im>imThresh);
    oneStdIm = median(im(indThresh));
    ratioIm(i) = oneStdIm/imMed;
end

plot(3.*(([1:numPts]-1)),ratioIm, 'r');
legend('45min10inten25au', 'nobluelight');
title('Pho4-NLSvar#11, Pho5-Venus');
xlabel('time (min)');
ylabel('ratio fluorescence intensity over bg (au)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd \\elsamad.ucsf.edu\data\instrumentation\microscope\SYC\20151016_Pho4Pho5_Dot6Arx1_NOlight\Dot6Arx1_BlueLight45min10inten25au\
% extract YFP files
files = dir('*YFP*');
nTimes = size(files,1);
% extract values
A=[];
for i = 1:length(files)
    [mat,tok,ext]=regexpi(files(i).name,'\d+','match', 'tokens');
    A(i) = str2num(mat{2});
end
[valA, indA] = sort(A);
filesO = files(indA); % ordered

ratioIm = [];
for i = 1:length(filesO);
    im = imread(filesO(i).name);
    im = double(im);
    imMed = median(im(:));
    imStd = std(im(:));
    imThresh = imMed+1.*imStd;
    indThresh = find(im>imThresh);
    oneStdIm = median(im(indThresh));
    ratioIm(i) = oneStdIm/imMed;
end

numPts = 30;
figure(3); plot(3.*(([1:numPts]-1)),ratioIm); hold on;

cd \\elsamad.ucsf.edu\data\instrumentation\microscope\SYC\20151016_Pho4Pho5_Dot6Arx1_NOlight\BlueLightNo\
% extract YFP files
files = dir('YFP_p2*');
nTimes = size(files,1);
% extract values
A=[];
for i = 1:length(files)
    [mat,tok,ext]=regexpi(files(i).name,'\d+','match', 'tokens');
    A(i) = str2num(mat{2});
end
[valA, indA] = sort(A);
filesO = files(indA); % ordered

ratioIm = [];
for i = 1:length(filesO);
    im = imread(filesO(i).name);
    im = double(im);
    imMed = median(im(:));
    imStd = std(im(:));
    imThresh = imMed+1.*imStd;
    indThresh = find(im>imThresh);
    oneStdIm = median(im(indThresh));
    ratioIm(i) = oneStdIm/imMed;
end

plot(3.*(([1:numPts]-1)),ratioIm, 'r');
legend('45min10inten25au', 'nobluelight');
title('Dot6-NLSvar#11, Arx1-Venus');
xlabel('time (min)');
ylabel('ratio fluorescence intensity over bg (au)');
%%
%% 20151015 YFP