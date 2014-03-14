JSO code for microscope control as of 12/12/12

-->designed for analysis of images collected by JSO microscope scripts

-->primary functions are
	(1)CovMapCell -- takes image series from microscope, finds cells in each image
	(2)AlignCellsIMCT -- takes cell list from CovMapCell, constructs full cell traces
		-also see aligncell, aligncellsIM, aligncellsKmeans (all functional but some tradeoffs, IMCT is probably the best although slowest overall)


		
Contact JSO if problems occur (JacobStewartOrnstein@gmail.com)

Git notes: repository ignores any image files (eg background images) to get those look at the orignal files on the FTP server under Resources/script
this is after v1.0