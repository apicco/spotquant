import numpy as np
from os import listdir as ls
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
	
	#images acquired with the CCD are on a 12 bit range and can be rescaled to 16 bit
	im_rescaled  =  rescale_intensity( im , in_range = ( 0 , 2**12-1 ) , out_range = 'uint16' )
	img_corrected , img_median = img_corr( im_rescaled , radius )

	return img_corrected , img_median

def dilation(input,iterations):
	iter=0
	while iter < iterations:
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
		print( 'Yen threshold: ' + str( threshold_value ) )
	elif ( threshold_value == None ) & ( algorithm == 'otsu' ) :
		threshold_value = filters.threshold_otsu( image )
		print( 'Otsu threshold: ' + str( threshold_value ) )
	
	image_threshold = np.zeros( shape = image.shape , dtype = image.dtype )
	image_threshold[ image > threshold_value ]=1

	image_eroded = erosion( image_threshold , 21 )
	
	return image_eroded

#--------------------------------------------------
#	Track the patches through the stack
#--------------------------------------------------
def measure_spot_intensities( image , patch_mask , ctrl_mask , cell_mask ):
	#Measure the spot intensities in image, identified by patch_mask, in the cell ctrl_mask, within or outside cell_mask

	#label the masks of the image to be quantified
	patch_label=label( patch_mask , connectivity=2 )
	ctrl_label=label( ctrl_mask , connectivity=2 )

	measurements = np.zeros( ctrl_label.max() )
	for i in range( ctrl_label.max() ):
		spot_at_the_edge = False #I need to define spot_at_the_edge here, because if match_boundary.size < 0 and spot_at_the_edge_of_movie = False, then the definition of spot_at_the_edge is missed
		for j in range( patch_label.max() ):
			spot_at_the_edge_of_movie = np.any([np.any(ctrl_label[0,:,:]==j+1),np.any(patch_label[patch_label.shape[0]-1,:,:]==j+1)]) #is the ctrl_mask in a frame at the beginnig or at the end of the movie?
			match_boundary = ctrl_label[ np.all( [ ctrl_label == i + 1 , cell_mask == 0 ] , axis = 0 ) ]
			if not ( match_boundary.size > 0 & spot_at_the_edge_of_movie ): #if the ctrl mask is not at the edge of the movie (spots can be at the bottom of the image) and the spot does not sit in the proximity of the cell edge, then blank the corresponding ctrl_mask so that the spot is not selected.
				ctrl_mask[ ctrl_label == i + 1 ] = 0
				measurements[ i ] = 0
			else :
				match = patch_label[ np.all( [ ctrl_label == i+1 , patch_label == j+1 ] , axis = 0 ) ]
				if match.size:
					#before proceeding with the quantification check that the 
					#spot to be quantified is not truncated at the beginning or 
					#at the end of the stack i.e.: that the label of the spot 
					#mask is not present either in any of the array entries of 
					#the first frame (0) or in any of those of the last frame 
					#(patch_label.size[0]-1). The presence of a spot at the edge
					#is recorded as a variable, but the spot is not excluded yet
					#to allow the algorith to recongize whether there are multiple
					#spots within the patch of the ctrl_mask. If the spot selection
					#would be done at this point one could not exculde that an other
					#GFP patch, which colocalize within the same patch of the 
					#ctrl_mask (hence too close), would not start or end on the
					#edges of the stack.
					spot_at_the_edge=np.any([np.any(patch_label[0,:,:]==j+1),np.any(patch_label[patch_label.shape[0]-1,:,:]==j+1)])
					#control it is not on the border of the frame
					if not spot_at_the_edge: 
						spot_at_the_edge=np.any([np.any(patch_label[:,0,:]==j+1),np.any(patch_label[:,patch_label.shape[1]-1,:]==j+1)])
					if not spot_at_the_edge: 
						spot_at_the_edge=np.any([np.any(patch_label[:,:,0]==j+1),np.any(patch_label[:,:,patch_label.shape[2]-1]==j+1)])
	
					#I want to avoid patches that are too close. Those patches
					#will have the same control_label, which is a dilated version of
					#patch_labels. If close patches
					#exist then the value measurements[i] would be assigned more 
					#than once. If this happens, at the second assignment the 
					#value of the measurement will be non-zero. Therefore the 
					#measurement[i] will be deleted (assigned 0) and the loop
					#over range(patch_labels.max()) will exit to proceed with the
					#next patch in the ctrl_mask.
					if not measurements[i]:
						frame_measurements=np.zeros( patch_label.shape[0] )
						for frame in range( patch_label.shape[0] ) :
							frame_measurements[ frame ] = image[ frame , : , : ][ patch_label[ frame , : , : ]  == j + 1 ].sum()
	
						#NB: frame_measurements.sum()  is equal to image[patch_label==j+1].sum()
						measurements[i]=frame_measurements.max()
					else:
						ctrl_mask[ ctrl_label == i + 1 ] = 0
						measurements[ i ] = 0
						break

		if spot_at_the_edge:
			ctrl_mask[ ctrl_label == i + 1 ] = 0
			measurements[ i ] = 0

	tiff.imsave( './tmp.tif' , ctrl_mask )
	
	return [ m for m in measurements if m > 0 ]

#--------------------------------------------------
#	Analyse the images
#--------------------------------------------------

def analysis(path_in,radius=17,file_pattern='GFP-FW',save_masks=True , only_membrane = False ):
	#define the array in which all the measurements will be stored
	output_measurements=np.zeros(0)
	images=ls(path_in)
	
	GFP_images=[img for img in images if file_pattern in img]
	
	for i in range(len(GFP_images)):
		GFP_im , GFP_median = load_image( path_in + GFP_images[i] , radius )
		print( path_in + GFP_images[i] )
		
		mask_patches = mask( GFP_im )
		ctrl_mask =  dilation( mask_patches , iterations = 2 )
		
		if only_membrane :
			cell_mask_tmp =  mask( GFP_median , algorithm = 'otsu' )
			cell_mask = morphology.erosion( cell_mask_tmp , morphology.ball(3) ).astype( cell_mask_tmp.dtype )
		else :
			cell_mask = ctrl_mask * 0
	
		#compute the mask of the image and save it
		patch_mask= dilation( mask_patches , iterations = 1 )
		output_measurements=np.concatenate((
			output_measurements,
			measure_spot_intensities(GFP_im , patch_mask , ctrl_mask , cell_mask )
			))
		
		#save the ctrl mask	
		if save_masks:

			if not os.path.exist( 'masks' ) :
				os.makedir( 'masks' )
			tiff.imsave( path_in + 'masks' + GFP_images[i].replace( file_pattern , '_CellMask.tif' ) , cell_mask )
			tiff.imsave( path_in + 'masks' + GFP_images[i].replace( file_pattern , '_CtrlMask.tif' ) , ctrl_mask )
			tiff.imsave( path_in + 'masks' + GFP_images[i].replace( file_pattern , '_PatchMask.tif' ) , patch_mask )

	return output_measurements

def experiment( path , target_name , reference_name = 'Nuf2' , target_median_radius = 6 , reference_median_radius = 17 , only_membrane = False , file_pattern = '.tif' ):
	
	reference = analysis( path + '/' + reference_name + '/' , radius = reference_median_radius , file_pattern = file_pattern , only_membrane  =  False )
	np.savetxt(path+'/'+reference_name+'_intensities.txt',reference)
	
	target = analysis( path + '/' + target_name + '/' , radius = target_median_radius , file_pattern = file_pattern , only_membrane  =  only_membrane )
	np.savetxt(path+'/'+target_name+'_intensities.txt',target)
	
	return(reference,target)
	

