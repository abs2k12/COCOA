# ==============================
#
# Thanks to https://github.com/ih64/XRB-phot/blob/master/xrbOpticalReduction.py
#
# ==============================
import itertools
import functools
import inspect
import logging
import random
import math
import numpy as np
import os
import subprocess
from subprocess import*
from param_rtt import *
from sim2obs_param import *
#from field_param import*
#import names
import glob, stat,linecache

def prepDAOPHOT(fwhm):
    with open('daophot.opt','w') as dao:

        dao.write('FWHM='+str(fwhm)+'\n')
        dao.write('FIT='+str(float(fwhm)*1.2)+'\n')
        dao.write('PSF='+str(float(fwhm)*2.7)+'\n')

        dao.write('READ=0.1\n')
        dao.write('GAIN='+str(gain)+'\n')

        dao.write('TH=5.0\n')

        dao.write('AN=1\n')

        dao.write('LOWBAD=15\n')

        dao.write('HIBAD='+str(satl)+'\n')

        dao.write('WATCH=-2.0\n')

        dao.write('VAR=0\n')

    with open('photo.opt','w') as photo:
        #set the photometry radius to the fwhm
        photo.write('A1='+str(float(fwhm)*1.5)+'\n')
#        photo.write('A2='+str(float(fwhm)*2.0)+'\n')
#        photo.write('A3='+str(float(fwhm)*3.0)+'\n')
#        photo.write('A4='+str(float(fwhm)*4.0)+'\n')
#        photo.write('A5='+str(float(fwhm)*5.0)+'\n')
#        photo.write('A6='+str(float(fwhm)*6.0)+'\n')
#        photo.write('A7='+str(float(fwhm)*7.0)+'\n')
#        photo.write('A8='+str(float(fwhm)*8.0)+'\n')
#        photo.write('A9='+str(float(fwhm)*9.0)+'\n')
#        photo.write('A9='+str(float(fwhm)*10.0)+'\n')
        #inner sky radius to 4*fwhm
        photo.write('IS='+str(float(fwhm)*4.0)+'\n')
        #outer sky radius to 11*fwhm
        photo.write('OS='+str(float(fwhm)*5.0)+'\n')

    with open('allstar.opt','w') as allstar:

        allstar.write('fit='+str(float(fwhm)*1.2)+'\n')

        allstar.write('isky='+str(float(fwhm)*2.0)+'\n')

        allstar.write('osky='+str(float(fwhm)*5.0)+'\n')

        allstar.write('watch=0\n')

        allstar.write('redet=1\n')

    return

def prepDAOPHOT2(fwhm):

    with open('daophot.opt','w') as dao:
        dao.write('FWHM='+str(fwhm)+'\n')
        dao.write('FIT='+str(float(fwhm)*2.0)+'\n')
        dao.write('PSF='+str(float(fwhm)*2.7)+'\n')

        dao.write('READ=0.1\n')
        dao.write('GAIN='+str(gain)+'\n')

        dao.write('TH=5.0\n')

        dao.write('AN=1\n')

        dao.write('LOWBAD=15\n')

        dao.write('HIBAD='+str(satl)+'\n')

        dao.write('WATCH=-2.0\n')

        dao.write('VAR=0\n')

    with open('photo.opt','w') as photo:
        #set the photometry radius to the fwhm
        photo.write('A1='+str(float(fwhm)*1.5)+'\n')
#        photo.write('A2='+str(float(fwhm)*2.0)+'\n')
#        photo.write('A3='+str(float(fwhm)*3.0)+'\n')
#        photo.write('A4='+str(float(fwhm)*4.0)+'\n')
#        photo.write('A5='+str(float(fwhm)*5.0)+'\n')
#        photo.write('A6='+str(float(fwhm)*6.0)+'\n')
#        photo.write('A7='+str(float(fwhm)*7.0)+'\n')
#        photo.write('A8='+str(float(fwhm)*8.0)+'\n')
#        photo.write('A9='+str(float(fwhm)*9.0)+'\n')
#        photo.write('A9='+str(float(fwhm)*10.0)+'\n')
        #inner sky radius to 4*fwhm
        photo.write('IS='+str(float(fwhm)*4.0)+'\n')
        #outer sky radius to 11*fwhm
        photo.write('OS='+str(float(fwhm)*5.0)+'\n')

    with open('allstar.opt','w') as allstar:
        allstar.write('fit='+str(float(fwhm)*2.0)+'\n')

        allstar.write('isky='+str(float(fwhm)*2.0)+'\n')

        allstar.write('osky='+str(float(fwhm)*5.0)+'\n')

        allstar.write('watch=0\n')

        allstar.write('redet=1\n')

    return

def prepDAOPHOT3(fwhm):

    with open('daophot.opt','w') as dao:
        dao.write('FWHM='+str(fwhm)+'\n')
        dao.write('FIT='+str(float(fwhm)*2.5)+'\n')
        dao.write('PSF='+str(float(fwhm)*2.7)+'\n')

        dao.write('READ=0.1\n')
        dao.write('GAIN='+str(gain)+'\n')
        dao.write('TH=5.0\n')
        dao.write('AN=1\n')
        dao.write('LOWBAD=15\n')
        dao.write('HIBAD='+str(satl)+'\n')
        dao.write('WATCH=-2.0\n')
        dao.write('VAR=0\n')

    with open('photo.opt','w') as photo:
        #set the photometry radius to the fwhm
        photo.write('A1='+str(float(fwhm)*1.5)+'\n')
#        photo.write('A2='+str(float(fwhm)*2.0)+'\n')
#        photo.write('A3='+str(float(fwhm)*3.0)+'\n')
#        photo.write('A4='+str(float(fwhm)*4.0)+'\n')
#        photo.write('A5='+str(float(fwhm)*5.0)+'\n')
#        photo.write('A6='+str(float(fwhm)*6.0)+'\n')
#        photo.write('A7='+str(float(fwhm)*7.0)+'\n')
#        photo.write('A8='+str(float(fwhm)*8.0)+'\n')
#        photo.write('A9='+str(float(fwhm)*9.0)+'\n')
#        photo.write('A9='+str(float(fwhm)*10.0)+'\n')
        #inner sky radius to 4*fwhm
        photo.write('IS='+str(float(fwhm)*4.0)+'\n')
        #outer sky radius to 11*fwhm
        photo.write('OS='+str(float(fwhm)*5.0)+'\n')

    with open('allstar.opt','w') as allstar:
        allstar.write('fit='+str(float(fwhm)*2.5)+'\n')
        allstar.write('isky='+str(float(fwhm)*2.0)+'\n')
        allstar.write('osky='+str(float(fwhm)*5.0)+'\n')
        allstar.write('watch=0\n')
        allstar.write('redet=1\n')

    return

def prepDAOPHOT4(fwhm):

    with open('daophot.opt','w') as dao:
        dao.write('FWHM='+str(fwhm)+'\n')
        dao.write('FIT='+str(float(fwhm)*2.5)+'\n')
        dao.write('PSF='+str(float(fwhm)*2.7)+'\n')

        dao.write('READ=0.1\n')
        dao.write('GAIN='+str(gain)+'\n')

        dao.write('TH=5.0\n')

        dao.write('AN=2\n')

        dao.write('LOWBAD=15\n')
        dao.write('HIBAD='+str(satl)+'\n')
        dao.write('WATCH=-2.0\n')
        dao.write('VAR=0\n')

    with open('photo.opt','w') as photo:
        #set the photometry radius to the fwhm
        photo.write('A1='+str(float(fwhm)*1.5)+'\n')
#        photo.write('A2='+str(float(fwhm)*2.0)+'\n')
#        photo.write('A3='+str(float(fwhm)*3.0)+'\n')
#        photo.write('A4='+str(float(fwhm)*4.0)+'\n')
#        photo.write('A5='+str(float(fwhm)*5.0)+'\n')
#        photo.write('A6='+str(float(fwhm)*6.0)+'\n')
#        photo.write('A7='+str(float(fwhm)*7.0)+'\n')
#        photo.write('A8='+str(float(fwhm)*8.0)+'\n')
#        photo.write('A9='+str(float(fwhm)*9.0)+'\n')
#        photo.write('A9='+str(float(fwhm)*10.0)+'\n')
        #inner sky radius to 4*fwhm
        photo.write('IS='+str(float(fwhm)*4.0)+'\n')
        #outer sky radius to 11*fwhm
        photo.write('OS='+str(float(fwhm)*5.0)+'\n')

    with open('allstar.opt','w') as allstar:
        allstar.write('fit='+str(float(fwhm)*2.5)+'\n')
        allstar.write('isky='+str(0.5)+'\n')
        allstar.write('osky='+str(float(fwhm)*5.0)+'\n')
        allstar.write('watch=0\n')
        allstar.write('redet=1\n')

    return

def daophotWrap(filename):
    with open('daophotGo.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename+'.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename+'.coo\n')
        f.write(filename+'.ap\n')
        f.write('pick\n')
        f.write(filename+'.ap\n')
        f.write('100,17.0\n')
        f.write(filename+'.list\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.list\n')
        f.write(filename+'.psf\n')

        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')

    #give executable permission to the shell script
    os.chmod('daophotGo.sh',0755)
    return

def daophot_cont(filename):
    with open('daophotGo2.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.lst2\n')
        f.write(filename+'.psf2\n')
        f.write('\n')
        f.write('\n')
        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')

    #give executable permission to the shell script
    os.chmod('daophotGo2.sh',0755)
    return

def daophot_cont2(filename):
    with open('daophotGo3.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.lst3\n')
        f.write(filename+'.psf2\n')
        f.write('\n')
        f.write('\n')
        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')

    #give executable permission to the shell script
    os.chmod('daophotGo3.sh',0755)
    return


def daophotWrap2(filename):
    '''creates a shell script that will call daophot to run find, pick, and psf
    on the image "filename"
    INPUT:
    filename: a string of the filename you want to run daophot on
    OUTPUT:
    daophotGo.sh: a shell script that wraps around daophot
    '''
    with open('daophotGo.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename+'.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename+'.coo\n')
        f.write(filename+'.ap\n')
        f.write('pick\n')
        f.write(filename+'.ap\n')
        #pick 20 stars, magnitude limit 20
        f.write('80,17.0\n')
        f.write(filename+'.list\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.list\n')
        f.write(filename+'.psf\n')

        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')

    #give executable permission to the shell script
    os.chmod('daophotGo.sh',0755)
    return

def daophotWrap3(filename):
    '''creates a shell script that will call daophot to run find, pick, and psf
    on the image "filename"
    INPUT:
    filename: a string of the filename you want to run daophot on
    OUTPUT:
    daophotGo.sh: a shell script that wraps around daophot
    '''
    with open('daophotGo.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename+'.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename+'.coo\n')
        f.write(filename+'.ap\n')
        f.write('pick\n')
        f.write(filename+'.ap\n')
        #pick 20 stars, magnitude limit 20
        f.write('60,17.0\n')
        f.write(filename+'.list\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.list\n')
        f.write(filename+'.psf\n')

        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')


    #give executable permission to the shell script
    os.chmod('daophotGo.sh',0755)
    return

def daophotWrap4(filename):
    '''creates a shell script that will call daophot to run find, pick, and psf
    on the image "filename"
    INPUT:
    filename: a string of the filename you want to run daophot on
    OUTPUT:
    daophotGo.sh: a shell script that wraps around daophot
    '''
    with open('daophotGo.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename+'.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename+'.coo\n')
        f.write(filename+'.ap\n')
        f.write('pick\n')
        f.write(filename+'.ap\n')
        #pick 20 stars, magnitude limit 20
        f.write('50,17.0\n')
        f.write(filename+'.lst\n')
        f.write('psf\n')
        f.write(filename+'.list\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.list\n')
        f.write(filename+'.psf\n')

        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')


    #give executable permission to the shell script
    os.chmod('daophotGo.sh',0755)
    return

def daophotWrap5(filename):
    '''creates a shell script that will call daophot to run find, pick, and psf
    on the image "filename"
    INPUT:
    filename: a string of the filename you want to run daophot on
    OUTPUT:
    daophotGo.sh: a shell script that wraps around daophot
    '''
    with open('daophotGo.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename+'.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename+'.coo\n')
        f.write(filename+'.ap\n')
        f.write('pick\n')
        f.write(filename+'.ap\n')
        #pick 20 stars, magnitude limit 20
        f.write('25,17.0\n')
        f.write(filename+'.list\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.list\n')
        f.write(filename+'.psf\n')

        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')


    #give executable permission to the shell script
    os.chmod('daophotGo.sh',0755)
    return

def daophotWrap6(filename):
    '''creates a shell script that will call daophot to run find, pick, and psf
    on the image "filename"
    INPUT:
    filename: a string of the filename you want to run daophot on
    OUTPUT:
    daophotGo.sh: a shell script that wraps around daophot
    '''
    with open('daophotGo.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename+'.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename+'.coo\n')
        f.write(filename+'.ap\n')
        f.write('pick\n')
        f.write(filename+'.ap\n')
        #pick 20 stars, magnitude limit 20
        f.write('10,17.0\n')
        f.write(filename+'.list\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.list\n')
        f.write(filename+'.psf\n')

        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')


    #give executable permission to the shell script
    os.chmod('daophotGo.sh',0755)
    return

def daophotWrap7(filename):
    '''creates a shell script that will call daophot to run find, pick, and psf
    on the image "filename"
    INPUT:
    filename: a string of the filename you want to run daophot on
    OUTPUT:
    daophotGo.sh: a shell script that wraps around daophot
    '''
    with open('daophotGo.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename+'.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename+'.coo\n')
        f.write(filename+'.ap\n')
        f.write('pick\n')
        f.write(filename+'.ap\n')
        #pick 20 stars, magnitude limit 20
        f.write('5,17.0\n')
        f.write(filename+'.list\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.list\n')
        f.write(filename+'.psf\n')

        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')

    #give executable permission to the shell script
    os.chmod('daophotGo.sh',0755)
    return

def daophotWrapcont1(filename):
    '''creates a shell script that will call daophot to run find, pick, and psf
    on the image "filename"
    INPUT:
    filename: a string of the filename you want to run daophot on
    OUTPUT:
    daophotGo.sh: a shell script that wraps around daophot
    '''
    with open('daophotGo.sh','w') as f:
        f.write('./daophot <<__DAOPHOT-END__\n')
        f.write('attatch '+filename+'\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename+'.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename+'.coo\n')
        f.write(filename+'.ap\n')
        f.write('pick\n')
        f.write(filename+'.ap\n')
        #pick 20 stars, magnitude limit 20
        f.write('100,17.0\n')
        f.write(filename+'.list\n')
        f.write('psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.list\n')
        f.write(filename+'.psf\n')

        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')

    #give executable permission to the shell script
    os.chmod('daophotGo.sh',0755)
    return

def allstarWrap(filename):
    '''creates a shells cript that will call allstar to do psf photometry on the image "filename"
    INPUT:
    filename: a string of the filename you want to do psfphotometry on
    OUTPUT:
    allstarGo.sh: a shell script that wraps around allstar
    '''
    with open('allstarGo.sh','w') as f:
        f.write('./allstar <<__ALLSTAR-END__\n')
        f.write('\n')
        f.write(filename+'\n')
        f.write(filename+'.psf2\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.als\n')
        f.write('\n')
        f.write('__ALLSTAR-END__\n')
        f.write('\n')

    #give executable permissions to shell script
    os.chmod('allstarGo.sh', 0777)
    pass

def allstarWrap2(filename):
    '''creates a shells cript that will call allstar to do psf photometry on the image "filename"
    INPUT:
    filename: a string of the filename you want to do psfphotometry on
    OUTPUT:
    allstarGo.sh: a shell script that wraps around allstar
    '''
    with open('allstarGo.sh','w') as f:
        f.write('./allstar <<__ALLSTAR-END__\n')
        f.write('\n')
        f.write(filename+'\n')
        f.write(filename+'.psf\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.als\n')
        f.write('\n')
        f.write('__ALLSTAR-END__\n')
        f.write('\n')

    #give executable permissions to shell script
    os.chmod('allstarGo.sh', 0777)
    pass

def allstarWrap3(filename):
    '''creates a shells cript that will call allstar to do psf photometry on the image "filename"
    INPUT:
    filename: a string of the filename you want to do psfphotometry on
    OUTPUT:
    allstarGo.sh: a shell script that wraps around allstar
    '''
    with open('allstarGo.sh','w') as f:
        f.write('./allstar <<__ALLSTAR-END__\n')
        f.write('\n')
        f.write(filename+'\n')
        f.write('0.0.fits.psf2\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.als\n')
        f.write('\n')
        f.write('__ALLSTAR-END__\n')
        f.write('\n')

    #give executable permissions to shell script
    os.chmod('allstarGo.sh', 0777)
    pass

def allstarWrap4(filename):
    '''creates a shells cript that will call allstar to do psf photometry on the image "filename"
    INPUT:
    filename: a string of the filename you want to do psfphotometry on
    OUTPUT:
    allstarGo.sh: a shell script that wraps around allstar
    '''
    with open('allstarGo.sh','w') as f:
        f.write('./allstar <<__ALLSTAR-END__\n')
        f.write('\n')
        f.write(filename+'\n')
        f.write('good.psf2\n')
        f.write(filename+'.ap\n')
        f.write(filename+'.als\n')
        f.write('\n')
        f.write('__ALLSTAR-END__\n')
        f.write('\n')

    #give executable permissions to shell script
    os.chmod('allstarGo.sh', 0777)
    pass


def remals1(filename):
    '''creates a shells cript that will call remals
    '''
    with open('remals.sh','w') as f:
        f.write('./remals << aaa'+'\n')

        f.write(filename+'.lst3'+'\n')
        f.write(filename+'.als\n')
        f.write(filename+'.rem\n')
        f.write('5\n')
        f.write('aaa')

    #give executable permissions to shell script
    os.chmod('remals.sh', 0777)
    pass

def remals2(filename):
    '''creates a shells cript that will call daophot/remals2
    '''
    with open('remals2.sh','w') as f:
        f.write('\rm -f '+filename+'.ap2'+'\n')
        f.write('\n')
        f.write('./daophot << aaa'+'\n')
        f.write('at '+filename+'\n')
        f.write('sub'+'\n')
        f.write(filename+'.psf2\n')
        f.write(filename+'.rem\n')
        f.write('n\n')
        f.write(filename+'s.fits\n')
        f.write('\n')
        f.write('at '+filename+'s.fits\n')
        f.write('\n')
        f.write('photo\n')
        f.write('photo2.opt\n')
        f.write('\n')
        f.write(filename+'.lst3'+'\n')
        f.write(filename+'.ap2'+'\n')
        f.write('\n')
        f.write('\n')
        f.write('\n')
        f.write('\n')
        f.write('sort'+'\n')
        f.write('4'+'\n')
        f.write(filename+'.ap2'+'\n')
        f.write(filename+'.srt'+'\n')
        f.write('\n')
        f.write('N\n')
        f.write('\n')
        f.write('exit\n')
        f.write('aaa\n')
        f.write('\n')
        f.write('mv -v '+filename+'.srt '+filename+'.ap2'+'\n')
        f.write('y\n')
        f.write('\n')
        f.write('\n')
        f.write('echo finished '+filename+'\n')

    #give executable permissions to shell script
    os.chmod('remals2.sh', 0777)
    pass

def phot2(fwhm):

    with open('photo2.opt','w') as photo:
        #set the photometry radius to the fwhm
        photo.write('A1='+str(float(fwhm)*1.5)+'\n')
        photo.write('A2='+str((float(fwhm)*1.5)+0.5)+'\n')
        photo.write('A3='+str((float(fwhm)*1.5)+1.0)+'\n')
        photo.write('A4='+str((float(fwhm)*1.5)+1.5)+'\n')
        photo.write('A5='+str((float(fwhm)*1.5)+2.75)+'\n')
        photo.write('A6='+str((float(fwhm)*1.5)+4.0)+'\n')
        photo.write('A7='+str((float(fwhm)*1.5)+5.25)+'\n')
        photo.write('A8='+str((float(fwhm)*1.5)+6.50)+'\n')
        photo.write('A9='+str((float(fwhm)*1.5)+7.80)+'\n')
        photo.write('AA='+str((float(fwhm)*1.5)+9.80)+'\n')
        photo.write('AB='+str((float(fwhm)*1.5)+11.5)+'\n')
        photo.write('AC='+str((float(fwhm)*1.5)+12.75)+'\n')
        #inner sky radius to 3*fwhm
        photo.write('IS='+str(float(fwhm)*3.0)+'\n')
        #outer sky radius to 4*fwhm
        photo.write('OS='+str(float(fwhm)*5.0)+'\n')
    return

def listair(filename):

    with open('list.air','w') as air:
        #set the photometry radius to the fwhm
        air.write(' '+filename+'.ap2                                1.001'+'\n')
    return

def listap(filename):

    with open('list.ap','w') as list:
        #set the photometry radius to the fwhm
        list.write(' '+filename+'.ap2'+'\n')
    return

def correction(filename):
    '''creates a shells cript that will call match2 and get aperture correction
    '''
    with open('correction.sh','w') as f:
        f.write('rm -f '+filename+'.cor'+'\n')
#        f.write('touch '+filename+'.cor'+'\n')
        f.write('\n')
        f.write('./match2 << aaa'+'\n')
        f.write('current.tot'+'\n')
        f.write(filename+'.als\n')
        f.write(filename+'.cor\n')
        f.write('10\n')
        f.write('y\n')
        f.write('y\n')
        f.write(filename+'.ref\n')
        f.write('aaa\n')
        f.write('\n')

        f.write('tail '+filename+'.cor >> current.cor')


    #give executable permissions to shell script
    os.chmod('correction.sh', 0777)
    pass
