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

The experiment folder must contain two subfolder. One subfoder contains the images of cells expressing the target protein tagged with your fluorophore of choice, the other contains the images of cells expressing the reference protein tagged with the same fluorophore. 


