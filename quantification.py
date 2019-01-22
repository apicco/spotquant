import measure_spot_intensities as msi
import numpy as np

images_path='.'
msi.experiment(path=images_path,target_name='Myo5_WT_24C',reference_name='Nuf2_WT_24C',GFP_pattern='100', target_median_radius = 6 , Nuf2_median_radius = 17 , only_membrane = False )
