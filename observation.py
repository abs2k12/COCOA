import itertools
import functools
import inspect
import logging
import random
import math
import numpy as np
import os
import subprocess
from param_rtt import *
#import names
from input_observation import *
import time

i=0
name_list = []


limitrarc = rtt*(180.0*3600.0/dist/math.pi)
#snapname = "projection0_10.dat"

LIMITING_RADIUS_X = limitrarc
LIMITING_RADIUS_Y = limitrarc

grid_count_y = int(LIMITING_RADIUS_Y / (valy * pixs)) + 1
grid_count_x = int(LIMITING_RADIUS_X / (valx * pixs)) + 1

bottom_left = -(valx * pixs * grid_count_x), -(valy * pixs * grid_count_y)
owd = os.getcwd()

os.mkdir("param-files")
os.mkdir("fits-files")
os.chdir("param-files")

for x in range(-grid_count_x + 1, grid_count_x):
    for y in range(-grid_count_y + 1, grid_count_y):
        print "Grid %s, %s" % (x, y)
        deco = y * (-pixs)*valy* OVERLAP_FACTOR
        raof = x * (-pixs)*valx * OVERLAP_FACTOR
        filename = "%s.%s" % (x, y)
        outfile = open(filename, 'w')
        outfile.write('NAXIS1     = '+ str(valx) + '       // image size'+'\n'+
        'NAXIS2     = '+ str(valy) + '       // image size' +'\n'+
        'DB-NX      = '+ str(posx) +'          // x-coordinate column in database' + '\n' +
        'DB-NY      = '+ str(posy) +'          // y-coordinate column in database'+ '\n' +
        'DB-NMAG    = '+ str(dmag) +'        // magnitude column in database' + '\n' +
        'OBJECT     = '+ str(obje) +'      // object name' + '\n' +
        'DISTANCE   = '+ str(dist) +'     // distance [parsecs]' +'\n'+
        'FILTER     = '+ str(filt) +''+'\n' +
        'PIXSCALE   = '+ str(pixs) +'      // pixel scale [arcseconds]' +'\n'+
        'GAIN       = '+ str(gain) +'      // [photons/ADU]'+'\n'+
        'SATLEVEL   = '+ str(satl) +'   // detector saturation level'+'\n'+
        'EXPOSURE   = '+ str(expo) +'        // fraction of a direct counts'+ '\n' +
        'SEEING     = '+ str(seei) +'        // [arcseconds]'+'\n' +
        'BACKGROUND = '+ str(back) +'       // background level' + '\n' +
        'PSF        = '+ str(psff) +'          // G-Gaussian, M-Moffat' + '\n' +
        'M_BETA     = '+ str(mbeta) +'       // beta parameter for Moffat' +'\n'+
        'NOISE      = '+ str(noise) +'          // 0: no noise, 1: Poisson noise'+'\n'+
        'RA_OFFSET  = '+ str(raof) +'        // RA offset from the center [arcseconds]' +'\n'+
        'DEC_OFFSET = '+ str(deco) + '        // DEC offset from the center [arcseconds]' +'\n'+
        'VERBOSE    = 2' +'\n' +
        'END')
        outfile.close()
        print ("deco: %s, raof: %s" % (deco, raof))
        outoff = open("offsets.dat", 'a')
        outoff.write(str(x)+ ' ' + str(y)+' '+ str(raof)+' ' + str(deco)+'\n')
        outoff.flush()

outoff.close()
time.sleep(2)
os.system("ls > names")
os.system('''sed '$d' names > allpar''')
os.system("rm names")
os.system("mv offsets.dat ../offsets.dat")
os.system("rm offsets.dat")

os.chdir(owd)
outfile2 = open("sim2obs_param.py", 'w')
outfile2.write("gain="+str(gain)+'\n'+"satl="+str(satl)+'\n'+ "dist="+ str(dist)+'\n'+ "pixs="+str(pixs)+'\n'+"sizx="+str(valx)+'\n'+"sizy="+str(valy)+'\n'+"overlap="+str(OVERLAP_FACTOR)+'\n')
outfile2.flush()
outfile2.close

outscript = open("auto-fits.sh", 'w')
outscript.write('#!/bin/bash'+'\n'+
'	i=1' +'\n'+
'        while IFS= read -r file' + '\n' +
'        do'+ '\n' +
'''                sim2obs/./sim2obs param-files/"$file" '''+ str(snapname)+''' fits-files/"$file".fits''' + '\n' +
'''done < param-files/"allpar"''')
outscript.flush()
outscript.close

os.system("bash auto-fits.sh")

#filename = "%s.%s" % (x, y)
# for names in open("nam1.dat", "r").readlines():
#     if names.endswith('\n'):
#     	names = names[:-1]
#     name_list.append(names)

# valp = 525
# valn = -525.0
# for names in name_list:
#     i = i+1
#     if i<=4:
#         outfile = open(str(names)+'.par', 'w')
#     	outfile.write(str(sizx)+'\n'+ str(sizy) +'\n'+ str(posx) + '\n' + str(posy) +'\n' + str(dmag) + '\n' + str(obje)
#                     + '\n' + str(dist) +'\n'+str(filt)+'\n' + str(pixs) +'\n'+str(gain)+'\n'+str(satl)+'\n'+str(expo)+ '\n' + str(seei)+'\n'
#                     + str(back) + '\n' + str(psff) + '\n' + str(mbet) +'\n'+str(nois)+'\n'
#                     + str(raof) +'\n'+ 'DEC_OFFSET = '+ str(valn) + '        // DEC offset from the center [arcseconds]' +'\n'
#                     + str(verb) +'\n' + str(end))
#     	outfile.flush()
#         valn=valn-525.0
#     else:
#     	outfile = open(str(names)+'.par', 'w')
#     	outfile.write(str(sizx)+'\n'+ str(sizy) +'\n'+ str(posx) + '\n' + str(posy) +'\n' + str(dmag) + '\n' + str(obje)
#                     + '\n' + str(dist) +'\n'+str(filt)+'\n' + str(pixs) +'\n'+str(gain)+'\n'+str(satl)+'\n'+str(expo)+ '\n' + str(seei)+'\n'
#                     + str(back) + '\n' + str(psff) + '\n' + str(mbet) +'\n'+str(nois)+'\n'
#                     + str(raof) +'\n'+ 'DEC_OFFSET = '+ str(valp) + '        // DEC offset from the center [arcseconds]' +'\n'
#                     + str(verb) +'\n' + str(end))
#     	outfile.flush()
#     	valp=valp+525.0
