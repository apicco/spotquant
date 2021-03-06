# Spotquant

Quantification of the fluorescence intensity in diffraction limited spots, imaged by fluorescence microscopy.

# Versions
1.04 - The threshold value used to compute the mask with the Yen algorithm is now computed on the image z max projection. That way the threshold values are smaller and the spot selection less stringent. Otherwise, dim spots were often discarded. 

# License
The software is distributed under the [GNU General Public License v3.0](https://github.com/apicco/spotquant/blob/master/LICENSE). However, we would appreciate if you cite [Picco, A., Kaksonen, M., _Precise tracking of the dynamics of multiple proteins in endocytic events_,  Methods in Cell Biology, Vol. 139, pages 51-68 (2017)](http://www.sciencedirect.com/science/article/pii/S0091679X16301546) if you use or modify this software. 

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
These `masks` folders contain all the masks used to identify and select the spots (\*\_PatchMask.tif and \*\_CtrlMask.tif), the images of the cells corrected for local background subtraction (\*\_BkgCorrected.tif); see Picco and Kaksonen, 2017 for details about the image background subtraction) as well as the median filtered image used to compute the background subtraction (\*\_Median.tif). These images are essential to asses the goodness of the parameters `target_median_radius` and `reference_median_radius`. A smaller median radius will leave the spots still visible on the \*\_Median.tif image. If that is the case, the radius must increase. Using the ImageJ preview, in `Process > Filters > Median...` is useful to find the ideal radius. A composite image of \*\_PatchMask.tif and its corresponding \*\_BkgCorrected.tif (in ImageJ: Color > Make Composite) will allow to identify the patches that have been selected:

![color_combine_example](https://github.com/apicco/spotquant/blob/master/example/example_of_patch_selection.png?raw=true)

(in red is the \*\_PatchMask.tif mask for that selected frame).

`save_masks` is an option to chose whether to save or not all the mask used to identify and measure spot intensities. It is set to `True` as default and creates a `masks` folder in the analised folder, which contains the images used for the analysis.

`measure_max_intensity_frame` defines how the fluorescence intensity of a spot is measured. If `True`, the fluorescence intensity is measured from the brightest frame in which each spot is imaged. That is the algorithm used in _Joglekar et al., 2006_ and _Picco et al., 2015_ and it is useful when the structures imaged are susceptible to move along the z-axis during the z-stack acquisition. With this option set as `True`, not all the photons captured by the objective are taken into account. When it is set to `False`, the fluorescence intensity is measured as the integral of the pixel values of each spot imaged in the different planes of the Z-stack.

# Analysis

`spotquant` has a basic function to estimate the number of molecules of the target protein spots, knowing how many molecules are present in the reference protein spots. 
Load the fluorescence intensities measured in the patches, which are saved as a .txt file, with

	import spotquant.measurespots as sq 
	foo_intensities = sq.load_quantification( "the_path_to_your_quantification/foo_intensities.txt" ) 
	foo_reference_intensities = sq.load_quantification( "the_path_to_your_quantification/foo_reference_intensities.txt" ) 

or use directly the outputs of `experiment`, which outputs both the reference and target protein intensities (the same that are saved into the .txt file).

The number of proteins can be estimated as

	foo_number_of_proteins = sp.quantify( foo_intensities , foo_reference_intensities , r_number ) 

where `r_number` is a tuple with two entries: the known number of proteins present in the reference patches, and its standard error. 

Refer to `example.py`, which analyses the images in the `example` folder, for a test using both `measure_max_intensity_frame = True` (default) and `measure_max_intensity_frame = False`. The number of target proteins will be (257.65, 7.21) and (277.86, 10.04).

# The function `quantify`

The fluorescence intensity values are computed as median (distributions are skewed). Their error is estimated with the MAD corrected for asymptotically normal consistency on the log transform of the raw fluorescence intensity values (used to aproximately conform to normality). 
The error associated with each fluorescence intensity value will thus be:

![error_MAD](https://latex.codecogs.com/gif.latex?%5Cbegin%7Balign*%7D%20l%20%26%3D%5Clog%28%20x%20%29%5Cquad%20x%3D%5Cleft%5C%7Bx_1%2C%5Cdots%2Cx_n%5Cright%5C%7D%5C%5C%20%5Chat%7Bl%7D%20%26%3D%5Ctext%7Bmedian%7D%28l%29%5C%5C%20%5Csigma_%7Bl%7D%26%3D%5Ctext%7BMAD%7D%28l%29/%5Csqrt%7Bn%7D%5C%5C%20%5Chat%7Bx%7D%20%26%3D%5Cexp%28%5Chat%7Bl%7D%29%3D%20%5Ctext%7Bmedian%7D%28%20x%20%29%5C%5C%20%5Csigma_%7B%5Chat%7Bx%7D%7D%20%26%3D%20%5Cexp%7B%28%5Chat%7Bl%7D%29%7D%5Csigma_%7Bl%7D%20%5Cend%7Balign*%7D)

where ![xx](https://latex.codecogs.com/gif.latex?x%3D%5Cleft%5C%7Bx_1%2C%5Cdots%2Cx_n%5Cright%5C%7D%5C%5C) are the measured fluorescence intensity values in each spot.

The error on the estimate of the number of molecules will be:

![average_molecules](https://latex.codecogs.com/gif.latex?e=n\frac{\hat{x}}{\hat{y}})

![error_molecules](https://latex.codecogs.com/gif.latex?\sigma_e=\sqrt{\left(\frac{\hat{x}}{\hat{y}}\sigma_n\right)^2+\left(n\frac{\hat{x}}{\hat{y}}\sigma_{\hat{x}}\right)^2+\left(n\frac{\hat{x}}{\hat{y}^2}\sigma_{\hat{y}}\right)^2})

where ![n](https://latex.codecogs.com/gif.latex?n) is the known number of molecules in the refrence spots, ![xx](https://latex.codecogs.com/gif.latex?\hat{x}) and ![yy](https://latex.codecogs.com/gif.latex?\hat{y}), are the median fluorescence intensity values measured for the reference and target protein respectively, and ![sn](https://latex.codecogs.com/gif.latex?\sigma_{\hat{x}})	, ![sr](https://latex.codecogs.com/gif.latex?\sigma_{\hat{y}}) are their respective error estimates derived from the log transformation described above. ![sr](https://latex.codecogs.com/gif.latex?\sigma_{n}) is the error estimate on the known number of molecules in the reference spot.

