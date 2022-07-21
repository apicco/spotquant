import spotquant.measurespots as sq
import pandas as pd

g , r = sq.spotquant( 'mdm34gfp-mdm12mcherry_imgs/' , file_pattern = 'tif' )

l = len( g )
if l == len( r ) :
    for i in range( l ) :
        if i == 0 :
            d = pd.DataFrame( [[ g[ i ] , r[ i ] ]] , columns = [ 'GFP' , 'RFP' ] )
        else :
            tmp = pd.DataFrame( [[ g[ i ] , r[ i ] ]] , columns = [ 'GFP' , 'RFP' ] )
            d = pd.concat( [ d , tmp ] )

d.to_csv( 'mdm34gfp_mdm12mcherry.csv' )
