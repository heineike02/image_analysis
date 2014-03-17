%USAGE
%
%This folder contains functions for analyzing some brightfield and
%fluorescence data acquired from the CellASIC microfluidic device. The
%scripts should be use to analyze cerevisiae cells.

%In the "ufluid_brightfield" folder, there are four .m files to note:
%   
%       imageanalysis_ufluid.m - is a script that runs the other three
%       functions.
%
%       loopthrough.m - is the function that loops through each frame of
%       the image series and stores the results of cell segmentation
%       obtained through bfanalyze.m in a structure
%
%       bfanalyze.m - is the function that actually does the image
%       segmentation and it simply takes a brightfield and a fluorescence
%       image at one time point and returns the relevant stats regarding
%       all the cells found in that time point
%
%       trackIDLu.m - is a function that takes the output from loopthrough
%       and feed it into the tracking algorithm, which is included in the
%       folder called "IDL_Particle_Tracking". Make sure to add this folder
%       to Matlab's search path
%
%Also in the folder is a folder "sampleplots", which contains some of the
%test plots I made in the process of doing the cell segmentation. The 
%graphs are for the first set of images, and they give you an idea of how 
%well Matlab tracked the cells for the first set of images.
