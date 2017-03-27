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

def true_anomaly(M,ecc,itermax=8):
    r"""
    Calculation of true and eccentric anomaly in Kepler orbits.
    
    ``M`` is the phase of the star, ``ecc`` is the eccentricity
    
    See p.39 of Hilditch, 'An Introduction To Close Binary Stars':
    
    Kepler's equation:
    
    .. math::
    
        E - e\sin E = \frac{2\pi}{P}(t-T)
        
    with :math:`E` the eccentric anomaly. The right hand size denotes the
    observed phase :math:`M`. This function returns the true anomaly, which is
    the position angle of the star in the orbit (:math:`\theta` in Hilditch'
    book). The relationship between the eccentric and true anomaly is as
    follows:
    
    .. math::
    
        \tan(\theta/2) = \sqrt{\frac{1+e}{1-e}} \tan(E/2)
    
    @parameter M: phase
    @type M: float
    @parameter ecc: eccentricity
    @type ecc: float
    @keyword itermax: maximum number of iterations
    @type itermax: integer
    @return: eccentric anomaly (E), true anomaly (theta)
    @rtype: float,float
    """
    # Initial value
    Fn = M + ecc*sin(M) + ecc**2/2.*sin(2*M)
    
    # Iterative solving of the transcendent Kepler's equation
    for i in range(itermax):
        F = Fn
        Mn = F-ecc*sin(F)
        Fn = F+(M-Mn)/(1.-ecc*cos(F))
        keep = F!=0 # take care of zerodivision
        if hasattr(F,'__iter__'):
            if np.all(abs((Fn-F)[keep]/F[keep])<0.00001):
                break
        elif (abs((Fn-F)/F)<0.00001):
            break
            
    # relationship between true anomaly (theta) and eccentric anomaly (Fn)
    true_an = 2.*arctan(sqrt((1.+ecc)/(1.-ecc))*tan(Fn/2.))
    
    return Fn,true_an  
