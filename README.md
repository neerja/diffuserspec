# diffuserspec

Diffuser Spectrometer Git Repo (MATLAB VERSION)

<img width="404" alt="Screen Shot 2024-10-24 at 4 22 07 PM" src="https://github.com/user-attachments/assets/5c598d17-bc27-4ec9-a2ad-5c942dc9b8bd">

Read the full paper at: https://opg.optica.org/ol/fulltext.cfm?uri=ol-48-2-323&id=524687

How to install this repository:

1. Clone the repo to your local computer. 
2. Download Raw Data folder from this link and place at diffuserspec/Raw Data: https://drive.google.com/drive/folders/1rjugHUgRvf3D8vrTu7i-7qbSZCdBH4ng?usp=sharing 

How to use this repository:

- Place raw data taken from experiments in DiffuserSpec/Raw Data/.  Be sure the filename has the dataset description AND date of acquisition (to help track).

- File size limit: X=0.025 GB (25 Mb).  Files larger than this MUST go in Raw Data and will be ignored during git commit.  Upload the large files to the Google Drive folder (above). 

- The hidden file .gitignore contains the large files that should be ignored during repo commit.  This file is HIDDEN by Finder and other system viewers.  But you can open via terminal. For more info: https://www.w3schools.com/git/git_ignore.asp

- Store any matFiles in Dataset matFiles.  This can include calibration files (if less than X GB).

- "Helper Functions" should contain any functions that will be used multiple times.

- Processing Code should contain bulk of code.  

- Results can contain png, fig, mp4 files of results. 

