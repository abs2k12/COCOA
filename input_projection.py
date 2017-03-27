# input file for projection
snapshot_name = "example_snap.dat"

seed = 10 #enter seed for random number generator, keeping the seed constant will produce the same initial projection


binary_velocities = 1

# Time input-intimes
#	-if control variable 'intimes' is set to 0 only one output file (projection0_*seed*.dat) is printed with projection of position and 
#velocities for single stars and binaries at t = 0
#	-if control variable 'intimes is set to 1 - multiple snapshots of binaries for 
#each time interval are printed, the final time and interval can be defined by the user in the variables below.
#   -if intimes = 2 then then multiple snapshots are also printed but the user can provide the times at which they are printe
# through an input file named 'time.dat'
# -output files produced when intimes is 1 or 2 are in the formation "projection0_*seed*_*time(in days)*.dat"
intimes = 0

#final time and interval can be specified in the units of days 
# can only when intimes is set to 1
finaltimes = 100.0
interval = 20.0

# Enable/Disable Rotation
# -if rotation = 0, output files give unprojected radius 
#to rotate the initial projection and multiple snapshots rotation needs to be equal to 1,
#roation angles alpha beta and gamma can be specified below
rotation = 0
#enter rotation angles  (all angles should be input in degrees) (can only be used if rotation = 1)     
alpha =  20.0	#rotation about x-axis
beta =   40.0	#rotation about y-axis
gamma =  60.0	 #rotation about z-axis

rc = 1e-1

#outputfiles produced when rotation is enabled are in the format "projection_*seed*_*alpha*_*beta*_*gamma*_time(in days)*.dat"


#----------------------------------------------------------------------------------------
#Calculate Surface Brightness Profile & Velocity Dispersion

sbpvdp = 0 #if sbpvdp is set to 1, Surface Brightness Profile and Velocity Dispersion is calculated

#parameters required for calculation of surface brightness profile and velocity dispersion

distance_to_sun= 5000.0 #the value of the source from the sun in the units of pc

reddening = 0.22 #value of reddening

maximum_magnitude = 2.285 #Maximum magnitude of the stars in the cluster


#----------------------------------------------------------------------------------------
#fit to king model

kingmodel = 0

option = 1              	# 0 - Trager and Gebhardt data, data from simulations SBP, VDP
distance_to_sun= 5000.0 	#the value of the source from the sun in the units of pc
av = 0.36           		# A_V  - redenning = 3*A_V
iline = 50             		# number of "good" lines in the input file (trager or model)
iline1 = 50             	# number of "good" lines in the input file (gephard)
input1 = 'sbpvdp.dat'		# input1 file - VDP (Gephard) or SBP and VDP
sigma1 = '5 30 7 1.3'     	# sigma (km/s), King scale radius r_c in arcsec, W0, m/l
sigma2 = '6 40 8 1.4'     	# sigma (km/s), King scale radius r_c in arcsec, W0, m/l

#----------------------------------------------------------------------------------------
#fit to king model_6 to make SBP

#kingsix = 0

#option_six = 1              	# 0 - Trager and Gebhardt data, data from simulations SBP, VDP
#distance_to_sun_six= 5000.0 	#the value of the source from the sun in the units of pc
#av_six = 0.36           		# A_V  - redenning = 3*A_V
#iline_six = 50             		# number of "good" lines in the input file (trager or model)
#sigma1_six = '0.2 18.0 6.0'     	#  x(1) = m/l, King scale radius r_c in arcsec, W0
#sigma2_six = '0.3 21.0 9.0'     	#  x(1) = m/l, King scale radius r_c in arcsec, W0



