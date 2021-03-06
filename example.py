import spotquant.measurespots as sq 

target , reference = sq.experiment( target_name = "Target_Protein" , path = "example/" , reference_name = 'Reference_Protein' , target_median_radius = 17 , reference_median_radius = 17 )

#-------------------------------------------------------#
#  If you already run sq.experiment, you can load its	#	
#  output uncommenting the two lines below:				#
#-------------------------------------------------------#
#target = sq.load_quantification( "example/Target_Protein_intensities.txt" )
#reference = sq.load_quantification( "example/Reference_Protein_intensities.txt" )

number_of_target_proteins = sq.quantify( target , reference , r_number = ( 80 , 0 ) )
print("The number of target proteins is ")
print( number_of_target_proteins )


target , reference = sq.experiment( target_name = "Target_Protein" , path = "example/" , reference_name = 'Reference_Protein' , target_median_radius = 17 , reference_median_radius = 17 , measure_max_intensity_frame = False , save_masks = False )
number_of_target_proteins = sq.quantify( target , reference , r_number = ( 80 , 0 ) )
print("The number of target proteins is ")
print( number_of_target_proteins )
