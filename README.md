# Spotquant

Quantification of the fluorescence intensity in diffraction limited spots, imaged by fluorescence microscopy.

# Installation

You need Python 3.X. To install the `spotquant` package it is sufficient to download it [here](https://github.com/apicco/spotquant/archive/master.zip).
Open the downloaded folder in a terminal and type

	sudo python setup.py install

The command `python` can be replaced by your Python 3.X installation of choice ( i.e. `python3.5`, `python3.7`, etc. etc.). 

# Use

It is sufficient to write a .py script with the following two lines

	import spotquant.measurespots as sq 
	sq.experiment( path  =  'a string containing the path to the experiment folder' , target_name = 'Foo1' , reference_name = 'Foo2' , target_median_radius  =  11 , reference_median_radius  =  17 , save_masks = True , measure_max_intensity_frame = True )

The experiment folder must contain two subfolders. 
One of these subfolders contains the images of cells expressing the target protein tagged with your fluorophore of choice, the other subfolder contains the images of cells expressing the reference protein tagged with the same fluorophore as the target proteins. 
It is important that the folders are named with the `target_name` and `reference_name` that are entered as variables in the `experiment` function. 
The script will output the fluorescence intensities as _.txt_ files in the folder identified by `path`. The script creates a `masks` folder in each of the two subfolders that contain the raw images. 
These `masks` folders contain all the masks used to identify and select the spots (\*\_PatchMask.tif and \*\_CtrlMask.tif) as well as the image of the cells corrected for local background subtraction (\*\_BkgCorrected.tif); see Picco and Kaksonen, 2017 for details about the image background subtraction) as well as the median filtered image used to compute the background subtraction (\*\_Median.tif). These images are essential to asses the goodness of the parameters `target_median_radius` and `reference_median_radius`. A smaller median radius will leave the spots still visible on the \*\_Median.tif image. If that is the case, the radius must increase. Using the ImageJ preview, in `Process > Filters > Median...` is useful to find the ideal radius.  

`save_mask` is an option to chose whether to save or not all the mask used to identify and measure spot intensities. It is set to `True` as default and creates a `masks` folder in the analised folder, which contains the images used for the analysis.

`measure_max_intensity_frame` defines how the fluorescence intensity of a spot is measured. If `True`, the fluorescence intensity is measured from the brightest frame in which each spot is imaged. That is the algorithm used in _Joglekar et al., 2006_ and _Picco et al., 2015_ and it is useful when the structures imaged are susceptible to move along the z-axis during the z-stack acquisition. With this option set as `True`, not all the photons captured by the objective are taken into account. When it is set to `False`, the fluorescence intensity is measured as the integral of the pixel values of each spot imaged in the different planes of the Z-stack.

