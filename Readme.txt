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

use BMH_data_analysis/BMH_20130716_analysis_klac.m as an example (has updated file structure).

Git notes: repository ignores any image files (eg background images) to get those look at the original files on the FTP server under Resources/script