BMH code for image analysis as of 17JUL2013
based on JSO code for image analysis as of 12/12/12

-->designed for analysis of images collected by JSO microscope scripts and micromanager

-->to see some examples of data analysis look in the BMH_data_analysis files

--> (1)CovMapCellsBMH -- takes image series from microscope, finds cells in each image
	(2)trackIDL_BMH replaces AlignCellsIMCT -- takes cell list from CovMapCell, constructs full cell traces
-also see aligncell, aligncellsIM, aligncellsKmeans (all functional but some tradeoffs, IMCT is probably the best although slowest overall)

Contact BMH if problems occur (heineike02@gmail.com)
Original version (1.0) created by JSO (JacobStewartOrnstein@gmail.com)

archive contains version 1.0 files that were not replaced in new version, but which were not necessary

use BMH_data_analysis/BMH_20140304_analysis_KL_vs_SC.m as an example (has updated file structure).

Git notes: repository ignores any image files (eg background images) to get those look at the original files on the FTP server under Resources/script

Desired features to add: 
Change image collectioncode output to sort in image order (01 v.s. 1 in filename), then adjust that within this structure. 

Use bright field to find cells and collect cell size
Spot cells with bright field
Collect Cell Size

For tracking adjust for final image shifts??

Filter out cells outside of a certain size range

make immune to os (use filesep) and computer (TrackIDL doesn't take ipdir argument)

Error handling for when no cells are detected.

To Do: 

Add information about circle files.  
storeim - is it necessary? 
- test std_thresh against nuclear localization

%Tracking parameters
%estimated maximum number of tracks

Are these even used? max_shift_rad
min_dist_thresh 

Test this vs new routine. %[maxx,maxy]=FindMaximaBMH(ima,px,py,close_max);


match_two_channels_Hung: 
%Make this function handle a situation where no cells are detected.

Noticed after implementing single channel / double channel functionality with cells as a channel (~7/30/2015) that a combinatorics error in the tracking began to 
occur.  

Test: run with and without changes and save xyzs from right before tracking.  

Why does load BG image take so long? 





