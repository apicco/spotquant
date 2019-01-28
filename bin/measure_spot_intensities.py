import numpy as np
from os import listdir as ls
from os import path
from os import mkdir
from skimage import filters #import filters
from skimage import morphology
from skimage.measure import label
from skimage.exposure import rescale_intensity
from skimage import img_as_float, img_as_uint
import skimage.external.tifffile as tiff
import skimage.io as io

#--------------------------------------------------
#	image correction
#--------------------------------------------------
def img_corr( img , radius ):

	#define imgage objects
	img_corrected  =  np.zeros( shape = img.shape , dtype = img.dtype )#where the image corrected will be stored
	img_median  =  np.zeros( shape = img.shape , dtype = img.dtype )#where the image median will be stored
	for i in range( img.shape[0] ):
	
		#compute the median of the cell for the frame i
		img_median[ i , : , : ] = filters.median( img[ i , : , : ] , morphology.disk( radius ) )
		
		#compute the image corrected without cytoplasmatic background. 
		#All maths are done on signed float dtype and converted in 'unsinged 16 bit' format
		img_corrected[ i , : , : ] = img_as_uint( 
			( img_as_float( img[ i , : , : ] ) - img_as_float( img_median[ i , : , : ] ) )
			)

	return img_corrected  ,  img_median

def load_image( path , radius = 17 ):

	im  =  tiff.imread( path )

	#images acquired with the CCD are on a 12 bit range and were rescaled to 16 bit. <- DEPRECATED:
	# rescale_intensity as been removed as it induces an error in the pixel value. Pixel values are 
	# integers but rescaling a 12-bit image into a 16-bit image requires some values to be float, which cannot be. 
	# Hence, rescaling is a crude approximation that should be avoided when quantifying fluorescence intensities.
	#im_rescaled  =  rescale_intensity( im , in_range = ( 0 , 2**12-1 ) , out_range = 'uint16' )

	img_corrected , img_median = img_corr( im , radius ) #subtract the local background

	return img_corrected , img_median 

def dilation( input , iterations ):

	iter=0

	while iter < iterations :

		input=morphology.dilation(input, morphology.ball(2)).astype(input.dtype)
		iter+=1 

	return input

def erosion( image , n ) :

	if n == 0:
		n = 1
		print("n set to 1; eroding pixels with no neighbor makes no sense")
	if n > 26:
		n = 26
		print("n set to 26; number of neighbor pixels cannot exceed 26")

	brush = np.array([
			[#0
				[
					[ 1 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#1
				[
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#2
				[
					[ 0 , 0 , 1 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#3
				[
					[ 0 , 0 , 0 ],
					[ 1 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#4
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#5
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 1 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#6
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 1 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#7
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#8
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 1 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#9
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 1 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#10
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 1 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#11
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 1 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#12
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 1 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#13
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 1 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#14
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 1 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#15
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 1 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#16
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 1 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#17
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 1 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#18
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#19
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 1 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#20
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 1 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#21
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#22
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 1 ],
					[ 0 , 0 , 0 ]] ,
				],
			[#23
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 1 , 0 , 0 ]] ,
				],
			[#24
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ]] ,
				],
			[#25
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 1 , 0 ],
					[ 0 , 0 , 0 ]] ,
				[
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 0 ],
					[ 0 , 0 , 1 ]] ,
				]
			])
					
	image_output = np.zeros( shape = image.shape , dtype = image.dtype )
	
	eroded_images=np.zeros(shape=image.shape,dtype=image.dtype)
	
	for j in range(26):
	
		tmp_image = morphology.binary_erosion( image , brush[j,:,:,:] )
		eroded_images = eroded_images + tmp_image
		
		image_output[ eroded_images > 26 - n ] = 1
	
	return image_output


def mask(image , threshold_value = None , algorithm = 'yen' ):
	#make a mask from the image using either Yen or Otsu thresholding
	
	if ( threshold_value == None ) & ( algorithm == 'yen' ) :

		threshold_value = filters.threshold_yen( image )
		print( '	Yen threshold: ' + str( threshold_value ) )

	elif ( threshold_value == None ) & ( algorithm == 'otsu' ) :
		
		threshold_value = filters.threshold_otsu( image )
		print( '	Otsu threshold: ' + str( threshold_value ) )
	
	image_threshold = np.zeros( shape = image.shape , dtype = image.dtype )
	image_threshold[ image > threshold_value ] = 1

	image_eroded = erosion( image_threshold , 21 )
	
	return image_eroded

#--------------------------------------------------
#	Track the patches through the stack
#--------------------------------------------------

def measure_spot_intensities( image , patch_mask , cell_mask ):

	#label the spot masks of the image
	patch_label = label( patch_mask , connectivity=2 )
	# ctrl_mask is used to asses if patches are 
	# too close to the image boundaries and 
	# exclude them. Its dilation of patch_mask by 1 ensures that 
	# no fluorescence intensity is left behind
	ctrl_mask =  dilation( patch_mask, iterations = 1 )
	ctrl_label =  label( ctrl_mask )

	measurements = np.array( [] , dtype = np.float64 )

	for i in range( patch_label.max() ):

		is_spot_at_the_edge = np.any( [
			np.any( patch_label[ 0 , : , : ] == i + 1 ) , #is the patch at the beginning of the stack
			np.any( patch_label[ : , 0 , : ] == i + 1 ) , #is the patch at one side of the stack
			np.any( patch_label[ : , : , 0 ] == i + 1 ) , #is the patch at one side of the stack
			np.any( patch_label[ patch_label.shape[ 0 ] - 1 , : , : ] == i + 1 ) , #is the patch at the end of the stack
			np.any( patch_label[ : , patch_label.shape[ 1 ] - 1 , : ] == i + 1 ) , #is the patch at the other side of the stack
			np.any( patch_label[ : , : , patch_label.shape[ 2 ] - 1 ] == i + 1 ) #is the patch at the other side of the stack
			] )
	

		if not is_spot_at_the_edge :

			# We want to exlude patches that are too close one to the other. Two patches are 
			# too close if their masks are 2 or less pixels one apart from the other. When the
			# mask of those patches is dilate by 1 iteration to create ctrl_mask, they regoins
			# will touch and they will be labelled as one regin in ctrl_label. To select whether
			# on patch has neighbor that are too close, we compute the mask of that patch alone
			# (called spot_mask) and we dilate. If the size of the dileted region that masks the
			# patch is equal to the size of its mask in ctrl_mask, then there are no neighbor
			# patches. A neighbor patch would infact increase the size of the mask in ctrl_mask.
			# As the masks are all 1. The sizes of the individual patch masks are computed with
			# sums.
			spot_mask  =  np.zeros( shape = patch_mask.shape , dtype = patch_mask.dtype )
			spot_mask[ patch_label == i + 1 ] = 1 
			spot_mask = dilation( spot_mask , iterations = 1 )

			ctrl_label_id = ctrl_label[ patch_label == i + 1 ][ 0 ] #store the label ID of the ctrl_mask that corresponds to the spot patch_label == i + 1
			is_the_patch_isolated = len( ctrl_mask[ ctrl_label == ctrl_label_id ] ) == len( spot_mask[ spot_mask == 1 ] )

			if is_the_patch_isolated :

				# Measure the average intensity of the patch. 
				measurements = np.append( measurements , np.mean( image[ patch_label == i + 1 ] ) )
			
#			else : 
#
#				patch_mask[ patch_label == i + 1 ] = 0	
#
#		else : 
#			
#			patch_mask[ patch_label == i + 1 ] = 0	
#
	return measurements , patch_mask , ctrl_mask , spot_mask

#--------------------------------------------------
#	Analyse the images
#--------------------------------------------------

def analysis(path_in , radius = 17 , file_pattern = 'GFP-FW' , save_masks = True , only_membrane = False ):

	#define the array in which all the measurements will be stored
	output_measurements = np.zeros(0)

	images = ls(path_in)
	
	GFP_images = [ img for img in images if file_pattern in img ]
	
	print( "Path: " + path_in )

	for i in range( len( GFP_images ) ) :

		# Load the image
		GFP_im , GFP_median = load_image( path_in + GFP_images[i] , radius )
		print( "Image: " + GFP_images[i] )
	
		# Compute a mask of the patches and of the cell.. 
		patch_mask = mask( GFP_im )

		cell_mask =  mask( GFP_median , algorithm = 'otsu' )

#		# In some cases, it is convenient to measure the fluorescence
#		# intensity of patches that are presente only on the membrante
#		if only_membrane :
#
#			cell_mask_tmp =  mask( GFP_median , algorithm = 'otsu' )
#			cell_mask = morphology.erosion( cell_mask_tmp , morphology.ball(3) ).astype( cell_mask_tmp.dtype )
#
#		else :
#			cell_mask = ctrl_mask * 0
#	
#		patch_mask = dilation( mask_patches , iterations = 1 )

		measurement , patch_mask , ctrl_mask , spot_mask = measure_spot_intensities(GFP_im , patch_mask , cell_mask )
		output_measurements=np.concatenate((
			output_measurements, measurement
			))
	
		#save the ctrl mask	
		if not path.exists( path_in + 'masks/' ) :
			mkdir( path_in + 'masks/' )

		#tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_CellMask.tif' ) , cell_mask )
		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_SpotMask.tif' ) , spot_mask )
		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_CtrlMask.tif' ) , ctrl_mask )
		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_PatchMask.tif' ) , patch_mask )
		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_GFPMedian.tif' ) , GFP_median )
		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_GFPBkgCorrected.tif' ) , GFP_im )

	return output_measurements

def experiment( path , target_name , reference_name = 'Nuf2' , target_median_radius = 6 , reference_median_radius = 17 , only_membrane = False , file_pattern = '.tif' ):
	
	reference = analysis( path + '/' + reference_name + '/' , radius = reference_median_radius , file_pattern = file_pattern , only_membrane  =  False )
	np.savetxt(path+'/'+reference_name+'_intensities.txt',reference)
	
	target = analysis( path + '/' + target_name + '/' , radius = target_median_radius , file_pattern = file_pattern , only_membrane  =  only_membrane )
	np.savetxt(path + '/' + target_name + '_intensities.txt',target)
	
	return(reference,target)
	

