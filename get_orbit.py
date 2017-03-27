# ==============================
#
# Thanks to Phoebe Project 2.0: http://phoebe-project.org/2.0a/docs/future/future/phoebe.dynamics.html
# ==============================

import itertools
import functools
import inspect
import logging
import numpy as np
from numpy import pi,sqrt,cos,sin,tan,arctan
from scipy.optimize import newton,bisect, fmin, brentq
from true_anomaly import true_anomaly


def get_orbit(times, period, ecc, sma, t0, per0=0., long_an=0., incl=0.,
              dpdt=0., deccdt=0., dperdt=0., mass_conservation=True,
              component='primary', t0type='periastron passage'):
    r"""
    Construct an orbit in observer coordinates.
    """
    # If t0 is not the time of periastron passage, convert it:
    if t0type == 'superior conjunction':
        phshift = 0.
        t0 = from_supconj_to_perpass(t0, period, per0, phshift=phshift)
    elif not t0type == 'periastron passage':
        raise ValueError('t0type needs to be one of "superior conjunction" or "periastron passage"')
    
    #-- if dpdt is non-zero, the period is actually an array, and the semi-
    #   major axis changes to match Kepler's third law (unless
    #   `mass_conservation` is set to False)
    if dpdt!=0:
        period_ = period
        period = dpdt*(times-t0) + period_
        if mass_conservation and not np.isscalar(period):
             sma = sma/period[0]**2*period**2
        elif mass_conservation:
             sma = sma/period_**2*period**2
    #-- if dperdt is non-zero, the argument of periastron is actually an
    #   array
    if dperdt!=0.:
        per0 = dperdt*(times-t0) + per0
    #-- if deccdt is non-zero, the eccentricity is actually an array
    if deccdt!=0:
        ecc = deccdt*(times-t0) + ecc
    #-- compute orbit
    n = 2*pi/period
    ma = n*(times-t0)
    E,theta = true_anomaly(ma,ecc)
    r = sma*(1-ecc*cos(E))
    PR = r*sin(theta)
    #-- compute rdot and thetadot
    l = r*(1+ecc*cos(theta))#-omega))
    L = 2*pi*sma**2/period*sqrt(1-ecc**2)
    rdot = L/l*ecc*sin(theta)#-omega)
    thetadot = L/r**2
    #-- the secondary is half an orbit further than the primary
    if 'sec' in component.lower():
        theta += pi
    #-- take care of Euler angles
    theta_ = theta+per0
    #-- convert to the right coordinate frame
    #-- create some shortcuts
    sin_theta_ = sin(theta_)
    cos_theta_ = cos(theta_)
    sin_longan = sin(long_an)
    cos_longan = cos(long_an)
    #-- spherical coordinates to cartesian coordinates. Note that we actually
    #   set incl=-incl (but of course it doesn't show up in the cosine). We
    #   do this to match the convention that superior conjunction (primary
    #   eclipsed) happens at periastron passage when per0=90 deg.
    x = r*(cos_longan*cos_theta_ - sin_longan*sin_theta_*cos(incl))
    y = r*(sin_longan*cos_theta_ + cos_longan*sin_theta_*cos(incl))
    z = r*(sin_theta_*sin(-incl))
    #-- spherical vectors to cartesian vectors, and then rotated for
    #   the Euler angles Omega and i.
    vx_ = cos_theta_*rdot - sin_theta_*r*thetadot
    vy_ = sin_theta_*rdot + cos_theta_*r*thetadot
    vx = cos_longan*vx_ - sin_longan*vy_*cos(incl) 
    vy = sin_longan*vx_ + cos_longan*vy_*cos(incl)
    vz = sin(-incl)*vy_
    #-- that's it!
    return (x,y,z),(vx,vy,vz),(theta_,long_an,incl)
