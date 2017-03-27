#Code to project binaries & single stars on the plane of the sky (with functionaility to print multiple snapshots for
#binaries and to rotate the projected distances and velocities)
#Functions get_orbit & true_anamoly taken from PyPhoebe 2.0

import itertools
import functools
import inspect
import logging
import random
import math
import numpy as np
import os
import subprocess
from numpy import pi,sqrt,cos,sin,tan,arctan
from scipy.optimize import newton,bisect, fmin, brentq
from get_orbit import get_orbit
from read_lines import read_lines
from input_projection import *
from constants import *
from subprocess import Popen
from subprocess import call
from string import *
import time

random.seed(seed)

if rotation == 1:
    a =  alpha*np.pi/180.
    b =   beta*np.pi/180.
    g =  gamma*np.pi/180.
    cosa = cos(a)
    sina = sin(a)
    cosb = cos(b)
    sinb = sin(b)
    cosg = cos(g)
    sing = sin(g)
else:
    alpha = 0
    beta = 0
    gamma = 0

rc = 1e-1

if intimes == 1:
    ftseconds = finaltimes * daytosec #converting observation time & interval in days to seconds
    intervalsec = interval * daytosec

def frange(x, y, jump): #defining frange to be used for control loop when intimes = 1
  while x < y:
    yield x
    x += jump

times = 0

infile = open(snapshot_name, "r") #The code requires snapshot.dat from MOCCA as input
data = infile.readlines() #reading data from lines of snapshot.dat
infile.close()

outfile1 = open('projection0_' + str(seed)+ '.dat', 'w') #output file in which projected snapshot for both single and binary stars will be printed
if rotation == 1:
    outfile2 = open('rot_projection0_' + str(seed)+ '.dat', 'w')
outfile_params = open('param_rtt.py','w')

def process(data):
    binfiles = {}
    rtt = 0.0
    for line_variables in read_lines(data): #reading data from snapshot.dat line by line in this iterative loop

        for var_name, variable in line_variables.iteritems():
            globals()[var_name] = variable


        rc0 = 0.0
        rc1 = rc
        irc = 100
        for j in range(1, 50):
            rc0 = (j-1)*rc
            rc1 = j*rc
            if r >= rc0 and r < rc1:
                 irc = j
                 break
        if rtt < r and r-rtt<=2.0:
          rtt = r



        costheta = random.uniform(-1., 1.)      #choosing costheta randomly between -1 and 1
        phi = random.uniform(0., 1.)*2*np.pi    #choosing the angle phi randomly between 0 and 2Pi

        sintheta = math.sqrt(1.-costheta**2)
        rp = r*sintheta                         #calculating projected distance
        cosphi = cos(phi)
        sinphi = sin(phi)
        rpx = rp*cosphi
        rpy = rp*sinphi
        rpz = r*costheta


        if rpz < 10**-15 and rpz > -10**-15: #to eliminate very small values generated from machine error
            rpz = 0.0
        if rpx < 10**-15 and rpx > -10**-15:
            rpx = 0.0
        if rpy < 10**-15 and rpy > -10**-15:
            rpy = 0.0

        eta = random.uniform(0., 1.)*2*np.pi   #choosing the angle eta randomly between 0 and 2Pi
        coseta = cos(eta)
        sineta = sin(eta)

        #projection of radial and tangential velocities for single stars and center of mass of binaries
        if vr > 0:
            vrp = -vr*costheta + vt*coseta*sintheta
            vtx1 = vr*sintheta + vt*coseta*costheta
            vty1 = vt*sineta
        else:
            vrp = vr*costheta + vt*coseta*sintheta
            vtx1 = -vr*sintheta + vt*coseta*costheta
            vty1 = vt*sineta
        vtx = vtx1*cosphi - vty1*sinphi
        vty = vtx1*sinphi + vty1*cosphi

        if vtx < 10**-15 and vtx > -10**-15:
            vtx = 0.0
        if vty < 10**-15 and vty > -10**-15:
            vty = 0.0
        if vrp < 10**-15 and vrp > -10**-15:
            vrp = 0.

        if rotation == 1: #if rotation is set to 1, then the projected x,y,z values for position and velocities are rotated using angles alpha,beta,gamma
            rrpx = (cosb*cosg)*rpx + (-cosb*sing)*rpy + sinb*rpz
            rrpy = (cosa*sing + sina*sinb*cosg)*rpx + (cosa*cosg - sina*sinb*sing)*rpy + (-sina*cosb)*rpz
            rrpz = (sina*sing-cosa*sinb*cosg)*rpx + (sina*cosg+cosa*sinb*sing)*rpy + (cosa*cosb)*rpz

            if rrpz < 10**-15 and rrpz > -10**-15:
               rrpz = 0.0
            if rrpx < 10**-15 and rrpx > -10**-15:
               rrpx = 0.0
            if rrpy < 10**-15 and rrpy > -10**-15:
               rrpy = 0.0

            rvtx = (cosb*cosg)*vtx + (-cosb*sing)*vty + sinb*vrp
            rvty = (cosa*sing + sina*sinb*cosg)*vtx + (cosa*cosg - sina*sinb*sing)*vty + (-sina*cosb)*vrp
            rvrp = (sina*sing-cosa*sinb*cosg)*vtx + (sina*cosg+cosa*sinb*sing)*vty + (cosa*cosb)*vrp

            if rvtx < 10**-15 and rvtx > -10**-15:
               rvtx = 0.0
            if rvty < 10**-15 and rvty > -10**-15:
               rvty = 0.0
            if rvrp < 10**-15 and rvrp > -10**-15:
               rvrp = 0.0


        if sm2 == 0:   #if line from snapshot.dat is a single star then projected distances and velocities are written in output.dat
            indb = 0   #indb is the first variable printed in outputfile it indicates whether star is single (indb=0) or in a binary (indb=1)
            periodday = 0
            t0day = 0
            per0 = 0
            long_an = 0
            incl = 0
            mb = mbv1 + mv1

            if rotation == 1: #if rotation is 1 printing a different output file in which rotated projected distances and velocities are given
                outfile2.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rrpx)  + ' ' + "{:16.10e}".format(rrpy)
                 + ' ' + "{:16.10e}".format(rrpz)
                 + ' ' + "{:16.10e}".format(rvrp) + ' ' + "{:16.10e}".format(rvtx) + ' ' + "{:16.10e}".format(rvty)
                 + ' '  + "{:7}".format(idd1)
                 + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                 + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                 + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                 + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv1)
                 + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                 + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                 + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc)+ ' ' + "{:16.10e}".format(costheta)
                 + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                 + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                 + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                outfile2.flush()

                outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx)  + ' ' + "{:16.10e}".format(rpy)
                 + ' ' + "{:16.10e}".format(rpz)
                 + ' ' + "{:16.10e}".format(vrp) + ' ' + "{:16.10e}".format(vtx) + ' ' + "{:16.10e}".format(vty)
                 + ' '  + "{:7}".format(idd1)
                 + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                 + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                 + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                 + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv1)
                 + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mb) + ' ' + "{:16.10e}".format(mbv1)
                 + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                 + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt)  + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                 + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                 + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                 + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                outfile1.flush()
            else:
                outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx)  + ' ' + "{:16.10e}".format(rpy)
                 + ' ' + "{:16.10e}".format(rpz)
                 + ' ' + "{:16.10e}".format(vrp) + ' ' + "{:16.10e}".format(vtx) + ' ' + "{:16.10e}".format(vty)
                 + ' '  + "{:7}".format(idd1)
                 + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                 + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                 + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                 + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv1)
                 + ' ' + "{:16.10e}".format(mbv1) + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mb) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                 + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                 + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                 + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                 + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                 + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                outfile1.flush()


        else: #if line from snapshot.dat is a binary star then the following procedure is carried out
#            ecrandom = = random.uniform(0., 1.)
#            ecc = ecrandom * ecrandom

            indb = 1
            mb1 = mbv1 + mv1
            mb2 = mbv2 + mv2
            totalmass = sm1+sm2  #calculating total mass
            aau = a*rsuntoau     #converting semimajor axis from snapshot.dat from solar radius to AU u
            periodyrs = math.sqrt(aau**3/totalmass) #calculating period of binary in days
            period = periodyrs*yearstosec #converting period to seconds to pass it to get_orbit properly
            amet = aau*automet #converting semimajor axis to meters to pass it to get_orbit
            a1 = sm2/totalmass*amet  #calculating semimajor axis of the components in the binary
            a2 = sm1/totalmass*amet
            if a1 > a2:           #to make suer that a1 is always the semimajor axis of the primary component
                sma1 = a2
                sma2 = a1
            else:
                sma1 = a1
                sma2 = a2

            t0 = -random.random()*period  #initial condition t0 (time of perisastron passage) randomly taken as a fraction of the period
            per0 = random.uniform(0., 1.)*np.pi     #argument of periastron (radian), randomly generated between 0 and Pi
            long_an = random.uniform(0., 1.)*np.pi  #longitude of ascending node (radian) randomly generated between 0 and Pi
            incl = random.uniform(0., 1.)*np.pi    #inclination angle (rad) randomly generated between 0 and Pi
            times = 0 #setting observation time to 0
            periodday = period*(1./daytosec)
            t0day = t0 * (1./daytosec)

            pos1,velo1,euler1 = get_orbit(times,period,ecc,sma1,t0,per0=per0, #calling get_orbit to get projected position and velocities from center of mass of binary
                          long_an=long_an,incl=incl,component='primary')
            pos2,velo2,euler2 = get_orbit(times,period,ecc,sma2,t0,per0=per0,
                         long_an=long_an,incl=incl,component='secondary')

            #converting units of positions (obtained from get_orbit) from meters to pc
            x1 = pos1[0] * mettopc
            y1 = pos1[1] * mettopc
            z1 = pos1[2] * mettopc

            x2 = pos2[0] * mettopc
            y2 = pos2[1] * mettopc
            z2 = pos2[2] * mettopc

            #translating positions to center of cluster cooridnates
            rpx1 = rpx + x1
            rpy1 = rpy + y1
            rpz1 = rpz + z1

            rpx2 = rpx + x2
            rpy2 = rpy + y2
            rpz2 = rpz + z2

            if rotation == 1: #rotating the position of the two stars (if rotation = 1)
                rrpx1 = (cosb*cosg)*rpx1 + (-cosb*sing)*rpy1 + sinb*rpz1
                rrpy1 = (cosa*sing + sina*sinb*cosg)*rpx1 + (cosa*cosg - sina*sinb*sing)*rpy1 + (-sina*cosb)*rpz1
                rrpz1 = (sina*sing-cosa*sinb*cosg)*rpx1 + (sina*cosg+cosa*sinb*sing)*rpy1 + (cosa*cosb)*rpz1

                rrpx2 = (cosb*cosg)*rpx2 + (-cosb*sing)*rpy2 + sinb*rpz2
                rrpy2 = (cosa*sing + sina*sinb*cosg)*rpx2 + (cosa*cosg - sina*sinb*sing)*rpy2 + (-sina*cosb)*rpz2
                rrpz2 = (sina*sing-cosa*sinb*cosg)*rpx2 + (sina*cosg+cosa*sinb*sing)*rpy2 + (cosa*cosb)*rpz2

                if rrpz1 < 10**-15 and rrpz1 > -10**-15:
                    rrpz1 = 0.0
                if rrpx1 < 10**-15 and rrpx1 > -10**-15:
                    rrpx1 = 0.0
                if rrpy1 < 10**-15 and rrpy1 > -10**-15:
                    rrpy1 = 0.0

                if rrpz2 < 10**-15 and rrpz2 > -10**-15:
                    rrpz1 = 0.0
                if rrpx2 < 10**-15 and rrpx2 > -10**-15:
                    rrpx2 = 0.0
                if rrpy2 < 10**-15 and rrpy2 > -10**-15:
                    rrpy2 = 0.0


            #converting units of velocities (obtained from get_orbit) from m/s to km/s
            vx1 = velo1[0] * mettokm
            vy1 = velo1[1] * mettokm
            vz1 = velo1[2] * mettokm

            vx2 = velo2[0] * mettokm
            vy2 = velo2[1] * mettokm
            vz2 = velo2[2] * mettokm

            #translating velocities to center of cluster cooridnates
            vrp1 = vrp + vz1
            vtx1 = vtx + vx1
            vty1 = vty + vy1

            vrp2 = vrp + vz2
            vtx2 = vtx + vx2
            vty2 = vty + vy2

            if rotation == 1: #rotating the velocities of the two stars (if rotation = 1)
                rvtx1 = (cosb*cosg)*vtx1 + (-cosb*sing)*vty1 + sinb*vrp1
                rvty1 = (cosa*sing + sina*sinb*cosg)*vtx1 + (cosa*cosg - sina*sinb*sing)*vty1 + (-sina*cosb)*vrp1
                rvrp1 = (sina*sing-cosa*sinb*cosg)*vtx1 + (sina*cosg+cosa*sinb*sing)*vty1 + (cosa*cosb)*vrp1

                if rvtx1 < 10**-15 and rvtx1 > -10**-15:
                    rvtx1 = 0.0
                if rvty1 < 10**-15 and rvty1 > -10**-15:
                    rvty1 = 0.0
                if rvrp1 < 10**-15 and rvrp1 > -10**-15:
                    rvrp1 = 0.0

                rvtx2 = (cosb*cosg)*vtx2 + (-cosb*sing)*vty2 + sinb*vrp2
                rvty2 = (cosa*sing + sina*sinb*cosg)*vtx2 + (cosa*cosg - sina*sinb*sing)*vty2 + (-sina*cosb)*vrp2
                rvrp2 = (sina*sing-cosa*sinb*cosg)*vtx2 + (sina*cosg+cosa*sinb*sing)*vty2 + (cosa*cosb)*vrp2

                if rvtx2 < 10**-15 and rvtx2 > -10**-15:
                    rvtx2 = 0.0
                if rvty2 < 10**-15 and rvty2 > -10**-15:
                    rvty2 = 0.0
                if rvrp2 < 10**-15 and rvrp2 > -10**-15:
                    rvrp2 = 0.0


                if binary_velocities == 0:

                    #rotating center of mass velocities
                    rcvtx = (cosb*cosg)*vtx + (-cosb*sing)*vty + sinb*vrp
                    rcvty = (cosa*sing + sina*sinb*cosg)*vtx + (cosa*cosg - sina*sinb*sing)*vty + (-sina*cosb)*vrp
                    rcvrp = (sina*sing-cosa*sinb*cosg)*vtx + (sina*cosg+cosa*sinb*sing)*vty + (cosa*cosb)*vrp

                    if rcvtx < 10**-15 and rcvtx > -10**-15:
                        rvtx = 0.0
                    if rcvty < 10**-15 and rcvty > -10**-15:
                        rcvty = 0.0
                    if rcvrp < 10**-15 and rcvrp > -10**-15:
                        rcvrp = 0.0

                    outfile2.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rrpx1)  + ' ' + "{:16.10e}".format(rrpy1)
                     + ' ' + "{:16.10e}".format(rrpz1)
                     + ' ' + "{:16.10e}".format(rcvrp) + ' ' + "{:16.10e}".format(rcvtx) + ' ' + "{:16.10e}".format(rcvty)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                     + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile2.flush()

                    outfile2.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rrpx2)  + ' ' + "{:16.10e}".format(rrpy2)
                     + ' ' + "{:16.10e}".format(rrpz2)
                     + ' ' + "{:16.10e}".format(rcvrp) + ' ' + "{:16.10e}".format(rcvtx) + ' ' + "{:16.10e}".format(rcvty)                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                     + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile2.flush()

                    outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx1)  + ' ' + "{:16.10e}".format(rpy1)
                     + ' ' + "{:16.10e}".format(rpz1)
                     + ' ' + "{:16.10e}".format(rcvrp) + ' ' + "{:16.10e}".format(rcvtx) + ' ' + "{:16.10e}".format(rcvty)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv1)
                     + ' ' + "{:16.10e}".format(mbv1) + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mb1) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile1.flush()
                    #writing projected data for star 2 in the binary to output.dat
                    outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx2)  + ' ' + "{:16.10e}".format(rpy2)
                     + ' ' + "{:16.10e}".format(rpz2)
                     + ' ' + "{:16.10e}".format(rcvrp) + ' ' + "{:16.10e}".format(rcvtx) + ' ' + "{:16.10e}".format(rcvty)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv2)
                     + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2) + ' ' + "{:16.10e}".format(mb2) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile1.flush()

                else:
                    outfile2.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rrpx1)  + ' ' + "{:16.10e}".format(rrpy1)
                     + ' ' + "{:16.10e}".format(rrpz1)
                     + ' ' + "{:16.10e}".format(rvrp1) + ' ' + "{:16.10e}".format(rvtx1) + ' ' + "{:16.10e}".format(rvty1)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                     + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile2.flush()

                    outfile2.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rrpx2)  + ' ' + "{:16.10e}".format(rrpy2)
                     + ' ' + "{:16.10e}".format(rrpz2)
                     + ' ' + "{:16.10e}".format(rvrp2) + ' ' + "{:16.10e}".format(rvtx2) + ' ' + "{:16.10e}".format(rvty2)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                     + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile2.flush()

                    outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx1)  + ' ' + "{:16.10e}".format(rpy1)
                     + ' ' + "{:16.10e}".format(rpz1)
                     + ' ' + "{:16.10e}".format(vrp1) + ' ' + "{:16.10e}".format(vtx1) + ' ' + "{:16.10e}".format(vty1)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv1)
                     + ' ' + "{:16.10e}".format(mbv1) + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mb1) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile1.flush()
                    #writing projected data for star 2 in the binary to output.dat
                    outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx2)  + ' ' + "{:16.10e}".format(rpy2)
                     + ' ' + "{:16.10e}".format(rpz2)
                     + ' ' + "{:16.10e}".format(vrp2) + ' ' + "{:16.10e}".format(vtx2) + ' ' + "{:16.10e}".format(vty2)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv2)
                     + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2) + ' ' + "{:16.10e}".format(mb2) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile1.flush()

            else:
                if binary_velocities == 0:
                    outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx1)  + ' ' + "{:16.10e}".format(rpy1)
                     + ' ' + "{:16.10e}".format(rpz1)
                     + ' ' + "{:16.10e}".format(vrp) + ' ' + "{:16.10e}".format(vtx) + ' ' + "{:16.10e}".format(vty)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv1)
                     + ' ' + "{:16.10e}".format(mbv1) + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mb1) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')

                    #writing projected data for star 2 in the binary to output.dat
                    outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx2)  + ' ' + "{:16.10e}".format(rpy2)
                     + ' ' + "{:16.10e}".format(rpz2)
                     + ' ' + "{:16.10e}".format(vrp) + ' ' + "{:16.10e}".format(vtx) + ' ' + "{:16.10e}".format(vty)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv2)
                     + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2) + ' ' + "{:16.10e}".format(mb2) +  ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile1.flush()
                else:
                    #writing projected data for star 1 in the binary to output.dat
                    outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx1)  + ' ' + "{:16.10e}".format(rpy1)
                     + ' ' + "{:16.10e}".format(rpz1)
                     + ' ' + "{:16.10e}".format(vrp1) + ' ' + "{:16.10e}".format(vtx1) + ' ' + "{:16.10e}".format(vty1)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv1)
                     + ' ' + "{:16.10e}".format(mbv1) + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mb1) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')

                    #writing projected data for star 2 in the binary to output.dat
                    outfile1.write(str(indb)+' ' + str(im)  + ' ' + "{:16.10e}".format(rpx2)  + ' ' + "{:16.10e}".format(rpy2)
                     + ' ' + "{:16.10e}".format(rpz2)
                     + ' ' + "{:16.10e}".format(vrp2) + ' ' + "{:16.10e}".format(vtx2) + ' ' + "{:16.10e}".format(vty2)
                     + ' '  + "{:7}".format(idd1)
                     + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                     + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                     + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                     + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv2)
                     + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2) + ' ' + "{:16.10e}".format(mb2) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                     + ' ' + "{:16.10e}".format(mi1) +  ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)
                     + ' ' + "{:16.10e}".format(r) + ' ' + "{:16.10e}".format(vr) + ' ' + "{:16.10e}".format(vt) + ' ' + "{:3}".format(irc) + ' ' + "{:16.10e}".format(costheta)
                     + ' ' + "{:16.10e}".format(phi) + ' ' + "{:16.10e}".format(eta)
                     + ' ' + "{:16.10e}".format(periodday) + ' '  + "{:16.10e}".format(t0day)
                     + ' ' + "{:16.10e}".format(per0) + ' ' + "{:16.10e}".format(long_an) + ' ' + "{:16.10e}".format(incl) + '\n')
                    outfile1.flush()

            if intimes == 1: #if control variable is set to 1 then following procedure is carried out
                i = 1
                filenumber = finaltimes/interval
                for times in frange(0,ftseconds, intervalsec): #observation time is controlled by this loop
                    if i <= filenumber:
                        timesdays = times*(1./daytosec)
                        if i not in binfiles:
                            if rotation == 1:
                                binfiles[i] = open('projection_'+ str(seed) + '_' + str(alpha) + '_' + str(beta) + '_' + str(gamma) + '_' + str(timesdays) + '.dat', 'a') #opening multiple output files to print projected data for binaries, depending on time interval
                                #writing header for binarypos files
                                header2 = "time="+str(timesdays)+ ' ' +"seed="+str(seed) + ' ' + "alpha="+str(alpha) + ' ' + "beta="+str(beta) + ' ' + "gamma="+str(gamma)
                                binfiles[i].write("# " + str(header2)+ '\n')
                                binfiles[i].write("{:11}".format("Time (days)") + ' ' + "{:9}".format("       im")  + ' ' + "{:16}".format("rpx")  + ' ' + "{:16.10}".format("rpy")
                                + ' ' + "{:16.10}".format("rpz")
                                + ' ' + "{:16.10}".format("vrp") + ' ' + "{:16.10}".format("vtx") + ' ' + "{:16.10}".format("vty")
                                + ' ' + "{:5}".format(" idd1")
                                + ' ' + "{:5}".format(" idd2") + ' ' +  "{:5}".format("  ikb") + ' ' + "{:5}".format("  ik1") + ' ' + "{:5}".format("  ik2") + ' ' + "{:16.10}".format("sm1")
                                + ' ' + "{:16.10}".format("sm2") + ' ' + "{:16.10}".format("slum1") + ' ' + "{:16.10}".format("slum2") + ' ' + "{:16.10}".format("rad1")
                                + ' ' + "{:16.10}".format("rad2") + ' ' + "{:16.10}".format("spin1") + ' ' + "{:16.10}".format("spin2") + ' ' + "{:7}".format("iinter1")
                                + ' ' + "{:7}".format("iinter2") + ' ' + "{:16.10}".format("a") + ' ' + "{:16.10}".format("ecc") + ' ' + "{:16.10}".format("mv")
                                + ' ' + "{:16.10}".format("mbv") + ' ' + "{:16.10}".format("mi") + ' ' + "{:16.10}".format("mv1") + ' ' + "{:16.10}".format("mbv1")
                                + ' ' + "{:16.10}".format("mi1") + ' ' + "{:16.10}".format("mv2") + ' ' + "{:16.10}".format("mbv2") + ' ' + "{:16.10}".format("mi2")+ '\n')
                            else:
                                binfiles[i] = open('projection0' + '_' + str(seed) + '_' + str(timesdays) + '.dat', 'a') #opening multiple output files to print projected data for binaries, depending on time interval
                                #writing header for binarypos files
                                header2 = "time="+str(timesdays)+ ' ' +"seed="+str(seed) + ' ' + "alpha="+str(alpha) + ' ' + "beta="+str(beta) + ' ' + "gamma="+str(gamma)
                                binfiles[i].write("# " + str(header2)+ '\n')
                                binfiles[i].write("{:11}".format("Time (days)") + ' ' + '    ' + "{:9}".format(" im")  + ' ' + "{:16}".format("rpx")  + ' ' + "{:16.10}".format("rpy")
                                + ' ' + "{:16.10}".format("rpz")
                                + ' ' + "{:16.10}".format("vrp") + ' ' + "{:16.10}".format("vtx") + ' ' + "{:16.10}".format("vty")
                                + ' ' + "{:5}".format(" idd1")
                                + ' ' + "{:5}".format(" idd2") + ' ' +  "{:5}".format("  ikb") + ' ' + "{:5}".format("  ik1") + ' ' + "{:5}".format("  ik2") + ' ' + "{:16.10}".format("sm1")
                                + ' ' + "{:16.10}".format("sm2") + ' ' + "{:16.10}".format("slum1") + ' ' + "{:16.10}".format("slum2") + ' ' + "{:16.10}".format("rad1")
                                + ' ' + "{:16.10}".format("rad2") + ' ' + "{:16.10}".format("spin1") + ' ' + "{:16.10}".format("spin2") + ' ' + "{:7}".format("iinter1")
                                + ' ' + "{:7}".format("iinter2") + ' ' + "{:16.10}".format("a") + ' ' + "{:16.10}".format("ecc") + ' ' + "{:16.10}".format("mv")
                                + ' ' + "{:16.10}".format("mbv") + ' ' + "{:16.10}".format("mi") + ' ' + "{:16.10}".format("mv1") + ' ' + "{:16.10}".format("mbv1")
                                + ' ' + "{:16.10}".format("mi1") + ' ' + "{:16.10}".format("mv2") + ' ' + "{:16.10}".format("mbv2") + ' ' + "{:16.10}".format("mi2")+ '\n')

                        pos1,velo1,euler1 = get_orbit(times,period,ecc,sma1,t0,per0=per0, #calling get_orbit to obtain positions & velocities
                                          long_an=long_an,incl=incl,component='primary')
                        pos2,velo2,euler2 = get_orbit(times,period,ecc,sma2,t0,per0=per0,
                                         long_an=long_an,incl=incl,component='secondary')

                        #converting units of positions (obtained from get_orbit) from meters to pc
                        x1 = pos1[0] * mettopc
                        y1 = pos1[1] * mettopc
                        z1 = pos1[2] * mettopc

                        x2 = pos2[0] * mettopc
                        y2 = pos2[1] * mettopc
                        z2 = pos2[2] * mettopc

                        #translating positions to center of cluster cooridnates

                        rpx1 = rpx + x1
                        rpy1 = rpy + y1
                        rpz1 = rpz + z2

                        rpx2 = rpx + x2
                        rpy2 = rpy + y2
                        rpz2 = rpz + z2

                        #converting units of velocities (obtained from get_orbit) from m/s to km/s
                        vx1 = velo1[0] * mettokm
                        vy1 = velo1[1] * mettokm
                        vz1 = velo1[2] * mettokm

                        vx2 = velo2[0] * mettokm
                        vy2 = velo2[1] * mettokm
                        vz2 = velo2[2] * mettokm

                        #translating velocities to center of cluster cooridnates

                        vrp1 = vrp + vz1
                        vtx1 = vtx + vx1
                        vty1 = vty + vy1

                        vrp2 = vrp + vz2
                        vtx2 = vtx + vx2
                        vty2 = vty + vy2

                        #converting units of period, semimajor axis etc for printout in binarypos*.dat
                        # timesdays = times*(1./daytosec)
                        # arsun1 = sma1 * (1./automet) * (1./rsuntoau)
                        # arsun2 = sma2 * (1./automet) * (1./rsuntoau)
                        # periodday = period*(1./daytosec)
                        # t0day = t0 * (1./daytosec)
                        if rotation == 1:
                            rrpx1 = (cosb*cosg)*rpx1 + (-cosb*sing)*rpy1 + sinb*rpz1
                            rrpy1 = (cosa*sing + sina*sinb*cosg)*rpx1 + (cosa*cosg - sina*sinb*sing)*rpy1 + (-sina*cosb)*rpz1
                            rrpz1 = (sina*sing-cosa*sinb*cosg)*rpx1 + (sina*cosg+cosa*sinb*sing)*rpy1 + (cosa*cosb)*rpz1

                            rrpx2 = (cosb*cosg)*rpx2 + (-cosb*sing)*rpy2 + sinb*rpz2
                            rrpy2 = (cosa*sing + sina*sinb*cosg)*rpx2 + (cosa*cosg - sina*sinb*sing)*rpy2 + (-sina*cosb)*rpz2
                            rrpz2 = (sina*sing-cosa*sinb*cosg)*rpx2 + (sina*cosg+cosa*sinb*sing)*rpy2 + (cosa*cosb)*rpz2

                            if rrpz1 < 10**-15 and rrpz1 > -10**-15:
                                rrpz1 = 0.0
                            if rrpx1 < 10**-15 and rrpx1 > -10**-15:
                                rrpx1 = 0.0
                            if rrpy1 < 10**-15 and rrpy1 > -10**-15:
                                rrpy1 = 0.0

                            if rrpz2 < 10**-15 and rrpz2 > -10**-15:
                                rrpz1 = 0.0
                            if rrpx2 < 10**-15 and rrpx2 > -10**-15:
                                rrpx2 = 0.0
                            if rrpy2 < 10**-15 and rrpy2 > -10**-15:
                                rrpy2 = 0.0

                            rvtx1 = (cosb*cosg)*vtx1 + (-cosb*sing)*vty1 + sinb*vrp1
                            rvty1 = (cosa*sing + sina*sinb*cosg)*vtx1 + (cosa*cosg - sina*sinb*sing)*vty1 + (-sina*cosb)*vrp1
                            rvrp1 = (sina*sing-cosa*sinb*cosg)*vtx1 + (sina*cosg+cosa*sinb*sing)*vty1 + (cosa*cosb)*vrp1


                            if rvtx1 < 10**-15 and rvtx1 > -10**-15:
                                rvtx1 = 0.0
                            if rvty1 < 10**-15 and rvty1 > -10**-15:
                                rvty1 = 0.0
                            if rvrp1 < 10**-15 and rvrp1 > -10**-15:
                                rvrp1 = 0.0

                            rvtx2 = (cosb*cosg)*vtx2 + (-cosb*sing)*vty2 + sinb*vrp2
                            rvty2 = (cosa*sing + sina*sinb*cosg)*vtx2 + (cosa*cosg - sina*sinb*sing)*vty2 + (-sina*cosb)*vrp2
                            rvrp2 = (sina*sing-cosa*sinb*cosg)*vtx2 + (sina*cosg+cosa*sinb*sing)*vty2 + (cosa*cosb)*vrp2


                            if rvtx2 < 10**-15 and rvtx2 > -10**-15:
                                rvtx2 = 0.0
                            if rvty2 < 10**-15 and rvty2 > -10**-15:
                                rvty2 = 0.0
                            if rvrp2 < 10**-15 and rvrp2 > -10**-15:
                                rvrp2 = 0.0

                            binfiles[i].write("{:11}".format(timesdays) + ' ' + "{:9}".format(im)  + ' ' + "{:16.10e}".format(rrpx1)  + ' ' + "{:16.10e}".format(rrpy1)
                                 + ' ' + "{:16.10e}".format(rrpz1)
                                 + ' ' + "{:16.10e}".format(rvrp1) + ' ' + "{:16.10e}".format(rvtx1) + ' ' + "{:16.10e}".format(rvty1)
                                 + ' ' + "{:7}".format(idd1)
                                 + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                                 + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                                 + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                                 + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                                 + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                                 + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)+ '\n')

                            binfiles[i].write("{:11}".format(timesdays) + ' ' + "{:9}".format(im)  + ' ' + "{:16.10e}".format(rrpx2)  + ' ' + "{:16.10e}".format(rrpy2)
                                 + ' ' + "{:16.10e}".format(rrpz2)
                                 + ' ' + "{:16.10e}".format(rvrp2) + ' ' + "{:16.10e}".format(rvtx2) + ' ' + "{:16.10e}".format(rvty2)
                                 + ' '  + "{:7}".format(idd1)
                                 + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                                 + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                                 + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                                 + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                                 + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                                 + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)+ '\n')
                            binfiles[i].flush()
                            i=i+1
                        else:
                            binfiles[i].write("{:11}".format(timesdays) + ' ' + "{:9}".format(im)  + ' ' + "{:16.10e}".format(rpx1)  + ' ' + "{:16.10e}".format(rpy1)
                             + ' ' + "{:16.10e}".format(rpz1)
                             + ' ' + "{:16.10e}".format(vrp1) + ' ' + "{:16.10e}".format(vtx1) + ' ' + "{:16.10e}".format(vty1)
                             + ' '  + "{:7}".format(idd1)
                             + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                             + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                             + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                             + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                             + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                             + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)+ '\n')

                            #printing projected data for star 2 in binarypos.dat
                            binfiles[i].write("{:11}".format(timesdays) + ' ' + "{:9}".format(im)  + ' ' + "{:16.10e}".format(rpx2)  + ' ' + "{:16.10e}".format(rpy2)
                             + ' ' + "{:16.10e}".format(rpz2)
                             + ' ' + "{:16.10e}".format(vrp2) + ' ' + "{:16.10e}".format(vtx2) + ' ' + "{:16.10e}".format(vty2)
                             + ' '  + "{:7}".format(idd1)
                             + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                             + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                             + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                             + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                             + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                             + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)+ '\n')
                            binfiles[i].flush()
                            i=i+1 #increment i which controls the file number for each time interval

            if intimes == 2: #if control variable is set to 1 then following procedure is carried out
                i = 1
                times_list = []
                for times in open("time.dat", "r").readlines():
                     times_list.append(float(times))
                filenumber = len(times_list)
                times_list = np.multiply(times_list, daytosec)
                #times_list = daytosec*np.array(times_list)
                for times in times_list: #observation time is controlled by this loop
                    if i <= filenumber:
                        timesdays = times*(1./daytosec)
                        if i not in binfiles:
                            if rotation == 1:
                                binfiles[i] = open('projection_' + str(seed) + '_' + str(alpha) + '_' + str(beta) + '_' + str(gamma) + '_' + str(timesdays) + '.dat', 'a') #opening multiple output files to print projected data for binaries, depending on time interval
                                #writing header for binarypos files
                                header2 = "time="+str(timesdays)+ ' ' +"seed="+str(seed) + ' ' + "alpha="+str(alpha) + ' ' + "beta="+str(beta) + ' ' + "gamma="+str(gamma)
                                binfiles[i].write("# " + str(header2)+ '\n')
                                binfiles[i].write("{:11}".format("Time (days)") + ' ' + '    ' + "{:9}".format(" im")  + ' ' + "{:16}".format("rpx")  + ' ' + "{:16.10}".format("rpy")
                                + ' ' + "{:16.10}".format("rpz")
                                + ' ' + "{:16.10}".format("vrp") + ' ' + "{:16.10}".format("vtx") + ' ' + "{:16.10}".format("vty")
                                + ' ' + "{:7}".format("idd1")
                                + ' ' + "{:7}".format("idd2") + ' ' +  "{:2}".format("ikb") + ' ' + "{:2}".format("ik1") + ' ' + "{:2}".format("ik2") + ' ' + "{:16.10}".format("sm1")
                                + ' ' + "{:16.10}".format("sm2") + ' ' + "{:16.10}".format("slum1") + ' ' + "{:16.10}".format("slum2") + ' ' + "{:16.10}".format("rad1")
                                + ' ' + "{:16.10}".format("rad2") + ' ' + "{:16.10}".format("spin1") + ' ' + "{:16.10}".format("spin2") + ' ' + "{:2}".format("iinter1")
                                + ' ' + "{:2}".format("iinter2") + ' ' + "{:16.10}".format("a") + ' ' + "{:16.10}".format("ecc") + ' ' + "{:16.10}".format("mv")
                                + ' ' + "{:16.10}".format("mbv") + ' ' + "{:16.10}".format("mi") + ' ' + "{:16.10}".format("mv1") + ' ' + "{:16.10}".format("mbv1")
                                + ' ' + "{:16.10}".format("mi1") + ' ' + "{:16.10}".format("mv2") + ' ' + "{:16.10}".format("mbv2") + ' ' + "{:16.10}".format("mi2")+ '\n')
                            else:
                                binfiles[i] = open('projection0_' + str(seed) + '_' + str(timesdays) + '.dat', 'a') #opening multiple output files to print projected data for binaries, depending on time interval
                                #writing header for binarypos files
                                header2 = "time="+str(timesdays)+ ' ' +"seed="+str(seed) + ' ' + "alpha="+str(alpha) + ' ' + "beta="+str(beta) + ' ' + "gamma="+str(gamma)
                                binfiles[i].write("# " + str(header2)+ '\n')
                                binfiles[i].write("{:11}".format("Time (days)") + ' ' + "{:9}".format("       im")  + ' ' + "{:16}".format("rpx")  + ' ' + "{:16.10}".format("rpy")
                                + ' ' + "{:16.10}".format("rpz")
                                + ' ' + "{:16.10}".format("vrp") + ' ' + "{:16.10}".format("vtx") + ' ' + "{:16.10}".format("vty")
                                + ' ' + "{:7}".format("idd1")
                                + ' ' + "{:7}".format("idd2") + ' ' +  "{:2}".format("ikb") + ' ' + "{:2}".format("ik1") + ' ' + "{:2}".format("ik2") + ' ' + "{:16.10}".format("sm1")
                                + ' ' + "{:16.10}".format("sm2") + ' ' + "{:16.10}".format("slum1") + ' ' + "{:16.10}".format("slum2") + ' ' + "{:16.10}".format("rad1")
                                + ' ' + "{:16.10}".format("rad2") + ' ' + "{:16.10}".format("spin1") + ' ' + "{:16.10}".format("spin2") + ' ' + "{:2}".format("iinter1")
                                + ' ' + "{:2}".format("iinter2") + ' ' + "{:16.10}".format("a") + ' ' + "{:16.10}".format("ecc") + ' ' + "{:16.10}".format("mv")
                                + ' ' + "{:16.10}".format("mbv") + ' ' + "{:16.10}".format("mi") + ' ' + "{:16.10}".format("mv1") + ' ' + "{:16.10}".format("mbv1")
                                + ' ' + "{:16.10}".format("mi1") + ' ' + "{:16.10}".format("mv2") + ' ' + "{:16.10}".format("mbv2") + ' ' + "{:16.10}".format("mi2")+ '\n')

                        pos1,velo1,euler1 = get_orbit(times,period,ecc,sma1,t0,per0=per0, #calling get_orbit to obtain positions & velocities
                                          long_an=long_an,incl=incl,component='primary')
                        pos2,velo2,euler2 = get_orbit(times,period,ecc,sma2,t0,per0=per0,
                                         long_an=long_an,incl=incl,component='secondary')

                        #converting units of positions (obtained from get_orbit) from meters to pc
                        x1 = pos1[0] * mettopc
                        y1 = pos1[1] * mettopc
                        z1 = pos1[2] * mettopc

                        x2 = pos2[0] * mettopc
                        y2 = pos2[1] * mettopc
                        z2 = pos2[2] * mettopc

                        #translating positions to center of cluster cooridnates

                        rpx1 = rpx + x1
                        rpy1 = rpy + y1
                        rpz1 = rpz + z2

                        rpx2 = rpx + x2
                        rpy2 = rpy + y2
                        rpz2 = rpz + z2

                        #converting units of velocities (obtained from get_orbit) from m/s to km/s
                        vx1 = velo1[0] * mettokm
                        vy1 = velo1[1] * mettokm
                        vz1 = velo1[2] * mettokm

                        vx2 = velo2[0] * mettokm
                        vy2 = velo2[1] * mettokm
                        vz2 = velo2[2] * mettokm

                        #translating velocities to center of cluster cooridnates

                        vrp1 = vrp + vz1
                        vtx1 = vtx + vx1
                        vty1 = vty + vy1

                        vrp2 = vrp + vz2
                        vtx2 = vtx + vx2
                        vty2 = vty + vy2

                        #converting units of period, semimajor axis etc for printout in binarypos*.dat
                        # timesdays = times*(1./daytosec)
                        # arsun1 = sma1 * (1./automet) * (1./rsuntoau)
                        # arsun2 = sma2 * (1./automet) * (1./rsuntoau)
                        # periodday = period*(1./daytosec)
                        # t0day = t0 * (1./daytosec)
                        if rotation == 1:
                            rrpx1 = (cosb*cosg)*rpx1 + (-cosb*sing)*rpy1 + sinb*rpz1
                            rrpy1 = (cosa*sing + sina*sinb*cosg)*rpx1 + (cosa*cosg - sina*sinb*sing)*rpy1 + (-sina*cosb)*rpz1
                            rrpz1 = (sina*sing-cosa*sinb*cosg)*rpx1 + (sina*cosg+cosa*sinb*sing)*rpy1 + (cosa*cosb)*rpz1

                            rrpx2 = (cosb*cosg)*rpx2 + (-cosb*sing)*rpy2 + sinb*rpz2
                            rrpy2 = (cosa*sing + sina*sinb*cosg)*rpx2 + (cosa*cosg - sina*sinb*sing)*rpy2 + (-sina*cosb)*rpz2
                            rrpz2 = (sina*sing-cosa*sinb*cosg)*rpx2 + (sina*cosg+cosa*sinb*sing)*rpy2 + (cosa*cosb)*rpz2

                            if rrpz1 < 10**-15 and rrpz1 > -10**-15:
                                rrpz1 = 0.0
                            if rrpx1 < 10**-15 and rrpx1 > -10**-15:
                                rrpx1 = 0.0
                            if rrpy1 < 10**-15 and rrpy1 > -10**-15:
                                rrpy1 = 0.0

                            if rrpz2 < 10**-15 and rrpz2 > -10**-15:
                                rrpz1 = 0.0
                            if rrpx2 < 10**-15 and rrpx2 > -10**-15:
                                rrpx2 = 0.0
                            if rrpy2 < 10**-15 and rrpy2 > -10**-15:
                                rrpy2 = 0.0

                            rvtx1 = (cosb*cosg)*vtx1 + (-cosb*sing)*vty1 + sinb*vrp1
                            rvty1 = (cosa*sing + sina*sinb*cosg)*vtx1 + (cosa*cosg - sina*sinb*sing)*vty1 + (-sina*cosb)*vrp1
                            rvrp1 = (sina*sing-cosa*sinb*cosg)*vtx1 + (sina*cosg+cosa*sinb*sing)*vty1 + (cosa*cosb)*vrp1


                            if rvtx1 < 10**-15 and rvtx1 > -10**-15:
                                rvtx1 = 0.0
                            if rvty1 < 10**-15 and rvty1 > -10**-15:
                                rvty1 = 0.0
                            if rvrp1 < 10**-15 and rvrp1 > -10**-15:
                                rvrp1 = 0.0

                            rvtx2 = (cosb*cosg)*vtx2 + (-cosb*sing)*vty2 + sinb*vrp2
                            rvty2 = (cosa*sing + sina*sinb*cosg)*vtx2 + (cosa*cosg - sina*sinb*sing)*vty2 + (-sina*cosb)*vrp2
                            rvrp2 = (sina*sing-cosa*sinb*cosg)*vtx2 + (sina*cosg+cosa*sinb*sing)*vty2 + (cosa*cosb)*vrp2


                            if rvtx2 < 10**-15 and rvtx2 > -10**-15:
                                rvtx2 = 0.0
                            if rvty2 < 10**-15 and rvty2 > -10**-15:
                                rvty2 = 0.0
                            if rvrp2 < 10**-15 and rvrp2 > -10**-15:
                                rvrp2 = 0.0

                            binfiles[i].write("{:11}".format(timesdays) + ' ' + "{:9}".format(im)  + ' ' + "{:16.10e}".format(rrpx1)  + ' ' + "{:16.10e}".format(rrpy1)
                                 + ' ' + "{:16.10e}".format(rrpz1)
                                 + ' ' + "{:16.10e}".format(rvrp1) + ' ' + "{:16.10e}".format(rvtx1) + ' ' + "{:16.10e}".format(rvty1)
                                 + ' ' + "{:7}".format(idd1)
                                 + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                                 + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                                 + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                                 + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                                 + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                                 + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)+ '\n')

                            binfiles[i].write("{:11}".format(timesdays) + ' ' + "{:9}".format(im)  + ' ' + "{:16.10e}".format(rrpx2)  + ' ' + "{:16.10e}".format(rrpy2)
                                 + ' ' + "{:16.10e}".format(rrpz2)
                                 + ' ' + "{:16.10e}".format(rvrp2) + ' ' + "{:16.10e}".format(rvtx2) + ' ' + "{:16.10e}".format(rvty2)
                                 + ' '  + "{:7}".format(idd1)
                                 + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                                 + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                                 + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                                 + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                                 + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                                 + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)+ '\n')
                            binfiles[i].flush()
                            i=i+1
                        else:
                            binfiles[i].write("{:11}".format(timesdays) + ' ' + "{:9}".format(im)  + ' ' + "{:16.10e}".format(rpx1)  + ' ' + "{:16.10e}".format(rpy1)
                             + ' ' + "{:16.10e}".format(rpz1)
                             + ' ' + "{:16.10e}".format(vrp1) + ' ' + "{:16.10e}".format(vtx1) + ' ' + "{:16.10e}".format(vty1)
                             + ' '  + "{:7}".format(idd1)
                             + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                             + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                             + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                             + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                             + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                             + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)+ '\n')

                            #printing projected data for star 2 in binarypos.dat
                            binfiles[i].write("{:11}".format(timesdays) + ' ' + "{:9}".format(im)  + ' ' + "{:16.10e}".format(rpx2)  + ' ' + "{:16.10e}".format(rpy2)
                             + ' ' + "{:16.10e}".format(rpz2)
                             + ' ' + "{:16.10e}".format(vrp2) + ' ' + "{:16.10e}".format(vtx2) + ' ' + "{:16.10e}".format(vty2)
                             + ' '  + "{:7}".format(idd1)
                             + ' ' + "{:7}".format(idd2) + ' ' +  "{:2}".format(ikb) + ' ' + "{:2}".format(ik1) + ' ' + "{:2}".format(ik2) + ' ' + "{:16.10e}".format(sm1)
                             + ' ' + "{:16.10e}".format(sm2) + ' ' + "{:16.10e}".format(slum1) + ' ' + "{:16.10e}".format(slum2) + ' ' + "{:16.10e}".format(rad1)
                             + ' ' + "{:16.10e}".format(rad2) + ' ' + "{:16.10e}".format(spin1) + ' ' + "{:16.10e}".format(spin2) + ' ' + "{:2}".format(iinter1)
                             + ' ' + "{:2}".format(iinter2) + ' ' + "{:16.10e}".format(a) + ' ' + "{:16.10e}".format(ecc) + ' ' + "{:16.10e}".format(mv)
                             + ' ' + "{:16.10e}".format(mbv) + ' ' + "{:16.10e}".format(mi) + ' ' + "{:16.10e}".format(mv1) + ' ' + "{:16.10e}".format(mbv1)
                             + ' ' + "{:16.10e}".format(mi1) + ' ' + "{:16.10e}".format(mv2) + ' ' + "{:16.10e}".format(mbv2) + ' ' + "{:16.10e}".format(mi2)+ '\n')
                            binfiles[i].flush()
                            i=i+1 #increment i which controls the file number for each time interval


    for binfile in binfiles.values(): #closing multiple binarypos*.dat
        binfile.close()
    outfile1.close() #closing output.dat
    if rotation == 1:
        outfile2.close() #closing output.dat

    if sbpvdp == 1:
        dtosun = str(distance_to_sun)
        dtosun.encode()
        red = str(reddening)
        red.encode()
        magmax = str(maximum_magnitude)
        magmax.encode()
        if rotation == 0:
            filename = str('projection0_' + str(seed)+ '.dat')
            filename.encode()
            fileoutname = str('sbpvdp0_'+ str(seed)+ '.dat')
            fileoutname.encode()
        if rotation == 1:
            filename = str('projection0_' + str(seed) + '_' + str(alpha) + '_' + str(beta) + '_' + str(gamma) + '_0.dat')
            filename.encode()
            fileoutname = str('sbpvdp_'+ str(seed) + '_' + str(alpha) + '_' + str(beta) + '_' + str(gamma) + '_0.dat')
            fileoutname.encode()

        p = Popen(['./sbp-vdp'], stderr=None, stdin=subprocess.PIPE, stdout = subprocess.PIPE)
        p.stdin.write(filename+'\n')
        p.stdin.flush()
        p.stdin.write(fileoutname+'\n')
        p.stdin.flush()
        p.stdin.write(dtosun + '\n')
        p.stdin.flush()
        p.stdin.write(red+'\n')
        p.stdin.flush()
        p.stdin.write(magmax+'\n')
        p.stdin.flush()

    time.sleep(5)

    if kingmodel == 1:
        awksbp = '''awk {print $4,$5}' sbpvdp0_'''+str(seed)+'''.dat > sbp.dat'''
        awkvdp = '''awk {print $4,$6}' sbpvdp0_'''+str(seed)+'''.dat > vdp.dat'''
        os.system(awksbp)
        os.system(awkvdp)
        opt = str(option)
        opt.encode()
        dtosun = str(distance_to_sun)
        dtosun.encode()
        avn = str(av)
        avn.encode()
        ilines = str(iline)
        ilines.encode()
        iline1s = str(iline1)
        iline1s.encode()
        sigma1s = str(sigma1)
        sigma1s.encode()
        sigma2s = str(sigma2)
        sigma2s.encode()
        input1s = str(input1)
        input1s.encode()
        outfile3 = open('fit-king5.conf', 'w')
        if rotation == 0:
            filename1 = str('sbp.dat')
            filename1.encode()
            filename2 = str('vdp.dat')
            filename2.encode()
            fileoutnamesbp = str('kingsvb0_'+ str(seed) + '_0.dat')
            fileoutnamesbp.encode()
            fileoutnamevdp = str('kingvdp0_'+ str(seed) + '_0.dat')
            fileoutnamevdp.encode()
            print opt,dtosun,av,iline,iline1,filename2,fileoutnamevdp,fileoutnamesbp,sigma1s,sigma2s,input1s
        if rotation == 1:
            awksbprot = '''awk {print $4,$5}' sbpvdp_'''+str(seed)+'''_'''+str(alpha)+'''_'''+str(beta)+'''_''' +str(gamma)+ '''_0.dat > sbp-rot.dat'''
            awkvdprot = '''awk {print $4,$6}' sbpvdp_'''+ str(seed) + '''_''' + str(alpha) + '''_''' + str(beta) + '''_''' + str(gamma) + '''_0.dat > vdp-rot.dat'''
            os.system(awksbprot)
            os.system(awkvdprot)
            filename1 = str('sbp-rot.dat')
            filename1.encode()
            filename2 = str('vdp-rot.dat')           
            filename2.encode()
            fileoutnamesbp = str('kingsvb0_'+ str(seed) + '_' + str(alpha) + '_' + str(beta) + '_' + str(gamma) + '_0.dat')
            fileoutnamesbp.encode()
            fileoutnamevdp = str('kingvdp0_'+ str(seed) + '_' + str(alpha) + '_' + str(beta) + '_' + str(gamma) + '_0.dat')
            fileoutnamevdp.encode()
        outfile3.write(opt+'\n'+ dtosun +'\n'+ avn + '\n' + ilines +'\n' + iline1s + '\n' + filename1 + '\n' + fileoutnamesbp
                    + '\n' + filename2 + '\n' + fileoutnamevdp +'\n' + sigma1s +'\n'+sigma2s)
        outfile3.flush()
        outfile3.close

        king5 = call(['./fit-king5'],stdin=subprocess.PIPE)

#    if kingsix == 1:
#        opt_six = str(option_six)
#        opt_six.encode()
#        dtosun_six = str(distance_to_sun_six)
#        dtosun_six.encode()
#        avn_six = str(av_six)
#        avn_six.encode()
#        ilines_six = str(iline_six)
#        ilines_six.encode()
#        sigma1s_six = str(sigma1_six)
#        sigma1s_six.encode()
#        sigma2s_six = str(sigma2_six)
#        sigma2s_six.encode()

#        outfile4 = open('fit-king6.conf', 'w')
#        if rotation == 0:
#            filename2 = str('sbpvdp0_'+ str(seed) + '.dat')
#            #filename2 = str('copysbpvdp.dat')
#            filename2.encode()
#            fileoutnamesbp = str('king_svb_'+ str(seed) + '_0.dat')
#            fileoutnamesbp.encode()

#        if rotation == 1:
#            filename2 = str('king_svb_'+ str(seed) + '.dat')
#            #filename2 = str('copysbpvdp.dat')
#            filename2.encode()
#            fileoutnamesbp = str('king_svb_'+ str(seed) + '_0.dat')
#            fileoutnamesbp.encode()
#        outfile4.write(opt_six+'\n'+ dtosun_six +'\n'+ avn_six + '\n' + ilines_six +'\n' + filename2 + '\n' + fileoutnamesbp
#                    + '\n' + sigma1s_six +'\n'+sigma2s_six)
#        outfile4.flush()
#        outfile4.close

#        king6 = call(['./fit-king6'],stdin=subprocess.PIPE)

# king6(data)
    limitr = math.sqrt(rtt**2/2.0)
    outfile_params.write("rtt=" + str(limitr) + '\n' + "snapname=" + '''"projection0_''' + str(seed)+ '''.dat"''')
    outfile_params.flush()
    outfile_params.close
process(data)


