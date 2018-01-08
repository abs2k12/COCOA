#Input Parameters for sim2obs to generate FITS images

valx = 2048       #image size
valy = 2048       #image size'
posx = 3          #x-coordinate column in database'
posy = 4          #y-coordinate column in database'
dmag = 26         #magnitude column in database'
obje = 'test'       #object name
dist = 5.0e3     #distance to cluster [parsecs]'
filt = 'V'
pixs = 0.25      #pixel scale [arcseconds]'
gain = 1.0       #[photons/ADU]'
satl = 99000.0    #detector saturation level'
expo = 50.0        #multiple of a direct counts'
seei = 0.8        #[arcseconds]'
back = 1.0       #background level'
psff = 'M'          # G-Gaussian, M-Moffat'
mbeta = 1.7        #beta parameter for Moffat'
noise = 1         # 0: no noise, 1: Poisson noise'
raof = 0.0        # RA offset from the center [arcseconds]'
deco = 0.0        # DEC offset from the center [arcseconds]'
OVERLAP_FACTOR = 0.97 #Defines the factor by which the offset values overlap when multiple fields are created for imaging the entire cluster. The value needs to be less than 1.0. A value of 0.95, would equal to a 5% overlap in adjoining images. 
