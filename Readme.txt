HES Lab code for image analysis as of 20150807
based on JSO code for image analysis as of 12/12/12

1. Designed for cell finding, nuclear intensity and localization analysis of images collected by JSO microscope scripts, micromanager, and metamorph. 

2. Nuclear localization based on ratio of top 5 pixels to median nuclear intensity (Ref. JSO BPAC paper, Elowitz?)

3. The code is kept in the image_processing.git repository on the HES lab server.  

4. Image Processing: 
	a. Processed image data for each sample (can be multiple fields of view) in an experiment is stored in a cell called all_tracks_vec, with a corresponding all_times_vec in which times for that well are stored.  posvec contains a cell array with the directory prefixes for all fields of view and image sets processed.  
	b. Run from gui: Open up the image_processing_GUI.m gui file and run it.  
	c. To run from a script use the image_processing_template.m from the BMH_data_analysis file. 

	Key functions: 
	(1)FindCells.m  -- finds individual cells in a particular image, calculates nuclear localization (nf) and median nuclear intensity (nmi). 
	(2)CovMapCells.m -- wrapper for analyzing data from a set of images from the same movie. 
	(3)trackIDL.m -- finds tracks for each cell.  
	(4)time_series_analysis.m -- wrapper to generate tracks for cells a set of movies from the same continuous experiment.  
	
Contact BMH if problems occur (heineike02@gmail.com)
Original version (1.0) created by JSO (JacobStewartOrnstein@gmail.com)

archive contains version 1.0 files that were not replaced in new version, but which were not necessary

Git notes: repository ignores any image files (eg background images) to get those look at the original files on the FTP server under Resources/script



To do list: 
- Create a test image and some testing code to validate findCells

- Find Cells: 
	- Troubleshoot overabundance of identified cells
	- why do we use std_thresh? 
		Remove duplicates in first findmaxima call
		Remove duplicates with the same nf and nmi? 
	- Test this vs new routine. %[maxx,maxy]=FindMaximaBMH(ima,px,py,close_max);
		- if necessary, filter after centering. 
	- What is storeim doing? 
	
- File name: times_from_umanager_metadata
	problem: assumes Z-stack is used for all micromanager datasets

- Alter JSO scope code to create output with numbers in numerical order.  Ensure that doesn't break this code. 

- Use bright field to find cells and collect cell size

- For tracking account for whole image shifts. 

- Filter out cells outside of a certain size range

- Error handling for when no cells are detected.
	%Make this function handle a situation where no cells are detected.
	
- documentation for data structure. 

- refine and provide better help for gui

- Add information about circle files.  

- Calibrate with nuclear localization images. 

- Are these even used? max_shift_rad, min_dist_thresh 

- Gui for image visualization? 







