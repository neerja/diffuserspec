README - DiffuserSpec figure generation code

Main folder: Projects\smartOCT\Publication\Articles\smartOCT_diffuseSpec\manuscript_Code_Data
	- contains all code and data for processing diffuserSpec data and generating figures 
	for diffuserSpec manuscript

TO RUN CODE: 1) open one of the FigureX_Main functions within [...\manuscript_Code_Data\Code] directory
	2) run code
	3) code will prompt you to select a working directory. Select the server location labeled 
		in the comments of the code (same as directory labeled in this readme under main folder. 
	4) code will run to completion 

Sub-folders:

	Code: contains 3 Main files for generating figure 2, 3 and 4 (labeled as FigureX_Main)
	
		subFunctions - contains subfunctions required for running the Main functions

	Data: main data files needed for Main files in Code folder

		SSTM_3D.tif - full 3D SSTM acquired through diffuserSpec calibration procedure (details in text)
		SSTM_background.tif - SSTM background for subtraction from SSTM_3D.tif

		broadbandData - data for broadband reconstructions
		calibrationFiles - precalibrated spectrometer data
		SSTM_gt - raw SSTM ground truth data (not used in Main files, just here if needed)

	