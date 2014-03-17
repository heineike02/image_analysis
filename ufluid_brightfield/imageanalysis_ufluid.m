%this script runs the image processing code

%NOTE: currently the image analysis script does well for the first set of
%images, but not as well (finds noncells) for the second set if using the
%same parameters. 
%to address this, simply change the options in the "imfindcircles"
%function in matlab. This option is currently hard coded in the "bfanalyze"
%script

%prep inputs for cell segmentation function
imdir = '/Users/susanychen/Documents/MATLAB/ufluid_brightfield/sample_images/';

images = struct(); files=dir(strcat(imdir,'*.tif'));
for i = 1:(length(files)./2);
    if length(files)./2 < 10;
        images(i).BF = strcat('glucgal2T0',num2str(i),'XY01C1.tif');
        images(i).FLUOR = strcat('glucgal2T0', num2str(i), 'XY01C2.tif');
    else
        images(i).BF = strcat('glucgal2T',num2str(i),'XY01C1.tif');
        images(i).FLUOR = strcat('glucgal2T', num2str(i), 'XY01C2.tif');
    end
end

%run and store cell segmentation
timecoursedata = loopthrough(imdir, images);

%prepare the input information for tracking function
%define inputs for track cells 
maxdisp = 4;

datainput = timecoursedata;

%time vector
times = [1:length(datainput)];

%add the path of the IDL particle tracking
addpath('/Users/susanychen/Documents/MATLAB/ufluid_brightfield/IDL_Particle_Tracking')

%track cells
tracks = trackIDLu(datainput,times, maxdisp);