from distutils.core import setup
from spotquant.measurespots import header

setup( name = 'spot quantification' ,
		version = str( header( printit = False ) ),
		description = 'Utilities to quantify spot intensities in images' ,
		author = 'Andrea Picco',
		author_email = 'andrea.picco@unige.ch',
		url = 'http://apicco.github.io/spotquant/',
		download_url = 'https://github.com/apicco/spotquant/archive/master.zip',
		packages = [ 'spotquant' ],
		license = 'The software is distributed under the terms of the GNU General Public License Version 3, June 2007. Trajalign is a free software and comes with ABSOLUTELY NO WARRANTY.'
		)
