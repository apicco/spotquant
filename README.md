# Spotquant

Quantification of the fluorescence intensity in diffraction limited spots, imaged by fluorescence microscopy.

# Installation

You need Python 3.X. To install the `spotquant` package it is sufficient to download it [here](https://github.com/apicco/spotquant/archive/master.zip).
Open the downloaded folder in a terminal and type

	sudo python setup.py install

python can be replaced by your Python 3.X installation of choice. 

# Use

It is sufficient to write a .py script with the following two lines

	import spotquant.measurespots as sq 
	sq.experiment( path  =  'a string containing the path to the experiment folder' , target_name = 'Mdm34' , reference_name = 'Cse4' , target_median_radius  =  6 , reference_median_radius  =  17 )

The experiment folder must contain two subfolder. One of those subfoders contains the images of cells expressing the target protein tagged with your fluorophore of choice, the other subfolder contains the images of cells expressing the reference protein tagged with the same fluorophore as the target proteins. It is important that the folders are named with the `target_name` and `reference_name` that are entered as variabiles in the `experiment` function. The script will output the fluorescence intensities as .txt files in the folder identified by `path`. The scrip creates a `masks` folder in each of the two subfolders that contain the raw images. These `masks` folders contain all the masks used to identify and select the spots (\*\_PatchMask.tif and \*\_CtrlMask.tif) as well as the image of the cells corrected for local backgourn subtraction (\*\_BkgCorrected.tif) ; see Picco and Kaksonen, 2017 for details about the image background subtraction) as well as the median filtered image used to compute the background subtraction (\*\_Median.tif). These images are essential to asses the goodness of the parameters `target_median_radius` and `reference_median_radius`. A smaller median radius will leave the spots still visible on the \*\_Median.tif image. If that is the case, the radius must increase. Using the ImageJ preview, in `Process > Filters > Median...` is useful to find the ideal radius.  


