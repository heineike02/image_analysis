function [circ, siz, rad, maxdisp, std_thresh] = species_magnification_params(species, op_amp, ipdir, maxdisp_1x)
%circ_params: gives parameters for input into time_series_analysis based on
%species and microscope optical amplification.  
%
%STD thresh:  Previously 0.2 worked well for SC images.  0.16 seemed to work best for k.lactis on 16JUL images
%May just want to change this to output circ filename.

if strcmp(op_amp,'1x')
    multiplier = 1.0;
    if strcmp(species,'KL')
        %load image of KLactis cell - this was made in BMH_20130916_analysis_klac.m
        load([ipdir, 'circKL_1x.mat'])
        circ_out = circ;
        siz = [15,15];
        rad = 4;
        std_thresh = 0.16;
    elseif strcmp(species,'SC')
        %load image of S.Cerevisiae cell (made by Jacob)
        load([ipdir,'circSC_1x_syc.mat'])
        circ_out = circ; 
        siz = [17,17];
        rad = 6;    
        std_thresh = 0.2;
    else 
        'Error - incorrect species name'
    end
    
elseif strcmp(op_amp,'1p5x')
    multiplier = 1.5;
    if strcmp(species,'KL')
        %load image of KLactis cell - this was made in BMH_20130916_analysis_klac.m
        load([ipdir, 'circKL_1p5x.mat'])
        circ_out = circ;
        siz = [18,18];
        rad = 8;
        std_thresh = 0.16;
    elseif strcmp(species,'SC')
        %load image of S.Cerevisiae cell (made by Jacob)
        load([ipdir,'circSC_1p5x_syc.mat'])
        circ_out = circ;
        siz = [25,25]; %[25,25];
        rad = 11;
        std_thresh = 0.2; %0.2;
    else
        'Error - incorrect species name'
    end
      
    
else
    'Error - incorrect optical amplification parameter'
end

maxdisp = maxdisp_1x*multiplier;

end

