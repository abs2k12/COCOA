import itertools
import functools
import inspect
import logging
import random
import math
import numpy as np
import os
import subprocess
from input_observation import seei
from subprocess import*
from param_rtt import *
from sim2obs_param import *
from read_off import read_lines
from photo_functions import *
#from field_param import*
#import names
import glob, stat,linecache

photo2 = str('photo2.opt')
photo2.encode()
list_air =str('list.air')
list_air.encode()
list_ap = str('list.ap')
list_ap.encode()
three = str('3')
three.encode()
growpar = str('0.9 0.0')
growpar.encode()
minmagdg = str('0.015')
minmagdg.encode()

owd = os.getcwd()
name_list = os.listdir("fits-files")


os.chdir(owd)
"""
os.chdir("fwhm_code")
os.system("./fwhm ../fits-files/0.0.fits > fwhm-stats")
pipe2 = Popen('''tail -n 1 fwhm-stats''', shell=True, stdout=PIPE, executable='/bin/bash').stdout
output2 = pipe2.readlines()
strip_list = [item.strip() for item in output2]
fwhm = [line.strip().split('  ')[1] for line in strip_list]
fwhm = float(fwhm[0])
"""
fwhm = seei / pixs


os.chdir(owd)
outfwhm = open("fwhm.dat","w")
outfwhm.write(str(fwhm)+' '+str(pixs)+' '+str(sizx)+' '+str(sizy)+' '+str(overlap))
outfwhm.flush()
outfwhm.close()

infile = open("offsets.dat", "r") #The code requires snapshot.dat from MOCCA as input
data = infile.readlines() #reading data from lines of snapshot.dat
infile.close()

os.system("cp daophot fits-files/")
os.system("cp allstar fits-files/")
os.chdir("fits-files")
os.system("cp ../fwhm_code/* .")
os.system("cp ../remals .")
os.system("cp ../match2 .")
os.system("cp ../daogrow2006 .")

def process(data):

    for line_variables in read_lines(data):

        for var_name, variable in line_variables.iteritems():
            globals()[var_name] = variable
        name=str(gridx)+"."+str(gridy)+".fits"
        print name
        image = name
        """
        os.system("./fwhm "+ str(name)+" > fwhm-stats")
        pipe2 = Popen('''tail -n 1 fwhm-stats''', shell=True, stdout=PIPE, executable='/bin/bash').stdout
        output2 = pipe2.readlines()
        strip_list = [item.strip() for item in output2]
        fwhm = [line.strip().split('  ')[1] for line in strip_list]
        fwhm = float(fwhm[0])
#        fwhm = 2.40
        """
        fwhm = seei / pixs
        print fwhm
        if fwhm < 0.0 or fwhm > 15.0:
          os.system("./fwhm 0.0.fits > fwhm-stats")
          pipe2 = Popen('''tail -n 1 fwhm-stats''', shell=True, stdout=PIPE, executable='/bin/bash').stdout
          output2 = pipe2.readlines()
          strip_list = [item.strip() for item in output2]
          fwhm = [line.strip().split('  ')[1] for line in strip_list]
          fwhm = float(fwhm[0])
#        os.system("touch 0.lst 0.nei")
        prepDAOPHOT(fwhm)

        daophotWrap(image)
        os.system("bash daophotGo.sh")
        if os.stat(str(image)+".psf").st_size == 0:
          os.system("rm ./"+str(image)+".coo")
          os.system("rm ./"+str(image)+".ap")
          os.system("rm ./"+str(image)+".list")
          os.system("rm ./"+str(image)+".psf")
          daophotWrap2(image)
          os.system("bash daophotGo.sh")
        if os.stat(str(image)+".psf").st_size == 0:
          os.system("rm ./"+str(image)+".coo")
          os.system("rm ./"+str(image)+".ap")
          os.system("rm ./"+str(image)+".list")
          os.system("rm ./"+str(image)+".psf")
          daophotWrap3(image)
          os.system("bash daophotGo.sh")
        if os.stat(str(image)+".psf").st_size == 0:
          os.system("rm ./"+str(image)+".coo")
          os.system("rm ./"+str(image)+".ap")
          os.system("rm ./"+str(image)+".list")
          os.system("rm ./"+str(image)+".psf")
          daophotWrap4(image)
          os.system("bash daophotGo.sh")
        if os.stat(str(image)+".psf").st_size == 0:
          os.system("rm ./"+str(image)+".coo")
          os.system("rm ./"+str(image)+".ap")
          os.system("rm ./"+str(image)+".list")
          os.system("rm ./"+str(image)+".psf")
          daophotWrap5(image)
          os.system("bash daophotGo.sh")
        if os.stat(str(image)+".psf").st_size == 0:
          os.system("rm ./"+str(image)+".coo")
          os.system("rm ./"+str(image)+".ap")
          os.system("rm ./"+str(image)+".list")
          os.system("rm ./"+str(image)+".psf")
          daophotWrap6(image)
          os.system("bash daophotGo.sh")
        if os.stat(str(image)+".psf").st_size == 0:
          os.system("rm ./"+str(image)+".coo")
          os.system("rm ./"+str(image)+".ap")
          os.system("rm ./"+str(image)+".list")
          os.system("rm ./"+str(image)+".psf")
          daophotWrap7(image)
          os.system("bash daophotGo.sh")
#        if os.stat(str(image)+".psf").st_size == 0:
#          os.system("rm ./"+str(image)+".coo")
#          os.system("rm ./"+str(image)+".ap")
#          os.system("rm ./"+str(image)+".list")
#          os.system("rm ./"+str(image)+".psf")
#          daophotWrapcont1(image)
#          os.system("bash daophotGo.sh")
#          os.system("cp exit.lst "+str(image)+".list")
#          os.system("cp exit.psf "+str(image)+".psf")
#          os.system("cp exit.nei "+str(image)+".nei")
#          os.system("rm exit.lst exit.psf exit.nei")

        if os.stat(str(image)+".psf").st_size == 0:
          os.system("rm ./"+str(image)+".coo")
          os.system("rm ./"+str(image)+".ap")
          os.system("rm ./"+str(image)+".list")
          os.system("rm ./"+str(image)+".psf")
          prepDAOPHOT2(fwhm)
          daophotWrap(image)
          os.system("bash daophotGo.sh")
        if os.stat(str(image)+".psf").st_size == 0:
          os.system("rm ./"+str(image)+".coo")
          os.system("rm ./"+str(image)+".ap")
          os.system("rm ./"+str(image)+".list")
          os.system("rm ./"+str(image)+".psf")
          prepDAOPHOT3(fwhm)
          daophotWrap(image)
          os.system("bash daophotGo.sh")
        os.system("mv ./*.lst ./"+str(image)+".lst2")
	os.system("mv ./*.nei ./"+str(image)+".nei2")
        daophot_cont(image)
        os.system("bash daophotGo2.sh")
        os.system("mv ./*.lst ./"+str(image)+".lst3")
	os.system("mv ./*.nei ./"+str(image)+".nei3")
        daophot_cont2(image)
        os.system("bash daophotGo3.sh")
        os.system("mv ./*.lst ./"+str(image)+".lst")
	os.system("mv ./*.nei ./"+str(image)+".nei")

        if os.stat(str(image)+".psf2").st_size == 0:
          allstarWrap2(image)
        else:
         allstarWrap(image)
#        if gridx==2 and gridy==1:
#          allstarWrap3(image)

#        allstarWrap(image)
        os.system("bash allstarGo.sh")
        awklim = '''awk 'NR > 3 { print }' '''+str(name)+'''.als'''
        pipe = Popen(awklim, shell=True, stdout=PIPE, executable='/bin/bash').stdout
        output = pipe.readlines()
        starnum=[]
        xposition=[]
        yposition=[]
        mag=[]
        magerror=[]
        midpoint = (sizx + sizy)/4.0

#        os.system("cp ?.lst "+image+".lst")
#        remals1(image)
#        os.system("bash remals.sh")
#        phot2(fwhm)
#        remals2(image)
#        os.system("cp ?s.fits current.fits")
#        os.system("bash remals2.sh")
#        listair(image)
#        listap(image)

#        p = Popen(['./daogrow2006'], stderr=None, stdin=subprocess.PIPE, stdout = subprocess.PIPE)#
#        p.stdin.write(photo2+'\n')
#        p.stdin.flush()
#        p.stdin.write('\n')
#        p.stdin.flush()
#        p.stdin.write(list_air+'\n')
#        p.stdin.flush()
#        p.stdin.write(list_ap+'\n')
#        p.stdin.flush()
#        p.stdin.write(three + '\n')
#        p.stdin.flush()
#        p.stdin.write(growpar+'\n')
#        p.stdin.flush()
#        p.stdin.write(minmagdg+'\n')
#        p.stdin.flush()

#        os.system("mv ./*.tot current.tot")

#        print "done-daogrow"
#        correction(image)
#        os.system("bash correction.sh")
#        os.system("tail -1 current.cor")
#        os.system("cp current.cor "+image+".cor")

        for row in output:
         Data = row.split()
         starnum.append(float(Data[0]))
         xposition.append(-offx+(((float(Data[1])-midpoint)*pixs)-pixs))
         yposition.append(-offy+(((float(Data[2])-midpoint)*pixs)-pixs))
         mag.append(float(Data[3]))
         magerror.append(float(Data[4]))
        with open('catalogue.dat','a') as f:
            lis=[xposition,yposition,mag,magerror,starnum]
            for x in zip(*lis):
                f.write("{0} {1} {2} {3} {4}\n".format(*x))
                f.flush()
        os.system("rm "+str(gridx)+"."+str(gridy)+"s.fits")
        os.system("rm ./"+str(image)+".coo")
        os.system("rm ./"+str(image)+".ap")
        os.system("rm ./"+str(image)+".lst")
        os.system("rm ./"+str(image)+".psf")
        os.system("rm ./*.lst")
        os.system("rm ./*.nei")
        os.system("rm ./*.lst2")
        os.system("rm ./*.nei2")
        os.system("rm current.fits")
        os.system("rm ?s.fits")
        os.system("rm ./*s.fits")
    f.close()
process(data)

os.system("cp ../clean_overlap .")
os.system("cp ../fwhm.dat .")
os.system("./clean_overlap")
os.system('''awk '!x[$0]++' overlap.dat > identical.dat''')
os.system("grep -F -v -f identical.dat all.dat > clean_catalogue.dat")
os.system("cp clean_catalogue.dat ../")
os.system("rm *.dat")

# for name in name_list:
#     image = name
#     prepDAOPHOT(fwhm)
#     daophotWrap(image)
#     allstarWrap(image)
#     os.system("bash daophotGo.sh")
#     os.system("bash allstarGo.sh")
#     awklim = '''awk 'NR > 3 { print }' '''+str(name)+'''.als'''
#     pipe = Popen(awklim, shell=True, stdout=PIPE, executable='/bin/bash').stdout
#     output = pipe.readlines()
#     starnum=[]
#     xposition=[]
#     yposition=[]
#     mag=[]
#     magerror=[]
#     for row in output:
#      Data = row.split()
#      starnum.append(float(Data[0]))
#      xposition.append(float(Data[1]))
#      yposition.append(float(Data[2]))
#      mag.append(float(Data[3]))
#      magerror.append(float(Data[4]))
#     with open('catalogue.dat','a') as f:
#         lis=[xposition,yposition,mag,magerror,starnum]
#         for x in zip(*lis):
#             f.write("{0} {1} {2} {3} {4}\n".format(*x))
#             f.flush()
#     os.system("rm *s.fits")

# f.close()
