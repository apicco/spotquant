import numpy as np
import warnings as wr
from skimage import filters #import filters
from skimage import morphology
from skimage import img_as_float, img_as_uint
import tifffile as tiff

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

	with wr.catch_warnings() :
		wr.simplefilter( "ignore" , category = UserWarning )
		img_corrected , img_median = img_corr( im , radius ) #subtract the local background

	return img_corrected , img_median 

def dilation( input , iterations ):

	iter=0

	while iter < iterations :

		input=morphology.dilation(input, morphology.ball(2)).astype(input.dtype)
		iter+=1 

	return input

def mask(image , threshold_value = None , algorithm = 'yen' ):
	#make a mask from the image using either Yen or Otsu thresholding
	
	if ( threshold_value == None ) & ( algorithm == 'yen' ) :

		threshold_value = filters.threshold_yen( np.max( image , axis = 0 ) )
		print( '	Yen threshold: ' + str( threshold_value ) )

	elif ( threshold_value == None ) & ( algorithm == 'otsu' ) :
		
		threshold_value = filters.threshold_otsu( np.max( image , axis = 0 ) )
		print( '	Otsu threshold: ' + str( threshold_value ) )
	
	image_threshold = np.zeros( shape = image.shape , dtype = image.dtype )
	image_threshold[ image > threshold_value ] = 1

	image_eroded = erosion( image_threshold , 21 )
	
	return image_eroded

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

