import numpy as np
from spotquant.imagehandling import img_corr, load_image, dilation , erosion , mask 
from os import listdir as ls
from os import path
from os import mkdir
from skimage.measure import label
import skimage.external.tifffile as tiff

def header( version = 1.00 , year = 2019 , printit = True ) :

	if printit :

		print('|-----------------------------------------------------|')
		print('| SpotQuant version ' + str( version ) + '\t\t\t\t      |')
		print('| Repository url: https://github.com/apicco/spotquant |')
		print('| Copyright: ' + str( year ) + ' Andrea Picco. \t\t      |')
		print('|-----------------------------------------------------|\n')
	
	elif not printit :

		return version

	else :

		raise AttributeError('Please, if you want to print the header (printit = True) or if you want to return the verion number only (printit = False).')

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
				measurements = np.append( measurements , np.sum( image[ patch_label == i + 1 ] ) )
			
			else : 

				patch_mask[ patch_label == i + 1 ] = 0	

		else : 
			
			patch_mask[ patch_label == i + 1 ] = 0	

	return measurements , patch_mask , ctrl_mask

#--------------------------------------------------
#	Analyse the images
#--------------------------------------------------

def spotquant(path_in , radius = 17 , file_pattern = 'GFP-FW' , save_masks = True ):

	header()

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

		measurement , patch_mask , ctrl_mask = measure_spot_intensities(GFP_im , patch_mask , cell_mask )
		output_measurements=np.concatenate((
			output_measurements, measurement
			))
	
		#save the ctrl mask	
		if not path.exists( path_in + 'masks/' ) :
			mkdir( path_in + 'masks/' )

		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_CtrlMask.tif' ) , ctrl_mask )
		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_PatchMask.tif' ) , patch_mask )
		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_GFPMedian.tif' ) , GFP_median )
		tiff.imsave( path_in + 'masks/' + GFP_images[i].replace( file_pattern , '_GFPBkgCorrected.tif' ) , GFP_im )

	return output_measurements

def experiment( path , target_name , reference_name = 'Nuf2' , target_median_radius = 6 , reference_median_radius = 17 , file_pattern = '.tif' ):
	
	reference = spotquant( path + '/' + reference_name + '/' , radius = reference_median_radius , file_pattern = file_pattern )
	np.savetxt(path+'/'+reference_name+'_intensities.txt',reference)
	
	target = spotquant( path + '/' + target_name + '/' , radius = target_median_radius , file_pattern = file_pattern )
	np.savetxt(path + '/' + target_name + '_intensities.txt',target)
	
	return( target , reference)
	

