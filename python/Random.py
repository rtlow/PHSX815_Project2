#! /usr/bin/env python

import numpy as np
import scipy.special as special

#################
# Random class
#################
# class that can generate random numbers
class Random:
    """A random number generator class"""

    # initialization method for Random class
    def __init__(self, seed = 5555):
        self.seed = seed
        self.m_v = np.uint64(4101842887655102017)
        self.m_w = np.uint64(1)
        self.m_u = np.uint64(1)
        
        self.m_u = np.uint64(self.seed) ^ self.m_v
        self.int64()
        self.m_v = self.m_u
        self.int64()
        self.m_w = self.m_v
        self.int64()

    # function returns a random 64 bit integer
    def int64(self):
        self.m_u = np.uint64(self.m_u * 2862933555777941757) + np.uint64(7046029254386353087)
        self.m_v ^= self.m_v >> np.uint64(17)
        self.m_v ^= self.m_v << np.uint64(31)
        self.m_v ^= self.m_v >> np.uint64(8)
        self.m_w = np.uint64(np.uint64(4294957665)*(self.m_w & np.uint64(0xffffffff))) + np.uint64((self.m_w >> np.uint64(32)))
        x = np.uint64(self.m_u ^ (self.m_u << np.uint64(21)))
        x ^= x >> np.uint64(35)
        x ^= x << np.uint64(4)
        with np.errstate(over='ignore'):
            return (x + self.m_v)^self.m_w

    # function returns a random floating point number between (0, 1) (uniform)
    def rand(self):
        return 5.42101086242752217E-20 * self.int64()

    # function returns a random integer (0 or 1) according to a Bernoulli distr.
    def Bernoulli(self, p=0.5):
        if p < 0. or p > 1.:
            return 1
  
        R = self.rand()

        if R < p:
            return 1
        else:
            return 0

    # function returns a random double (0 to infty) according to an exponential distribution
    def Exponential(self, beta=1.):
      # make sure beta is consistent with an exponential
      if beta <= 0.:
        beta = 1.

      R = self.rand();

      while R <= 0.:
        R = self.rand()

      X = -np.log(R)/beta

      return X
      
      
    # function returns a random double according to a normal distribution
    # uses the Ratio-of-Uniforms method described in Section 7.3.8
    # and implementation from Section 7.3.9
    def Normal(self, mu=1., sig=1.):
        
        # trick to emulating a do-while construct in Python

        u = 0.
        v = 0.
        x = 0.
        y = 0.
        q = 0.
        while True:

            u = self.rand()
            v = 1.7156*(self.rand() - 0.5)
            x = u - 0.449871
            y = np.abs(v) + 0.386595
            q = np.square(x) + y*(0.19600*y-0.25472*x)

            if not ((q > 0.27597) and ((q > 0.27846) or (np.square(v) > -4.*np.log(u)*np.square(u)))):
                break

        return mu + sig*v/u


    # function returns a random integer according to a Poisson distribution
    # See section 7.3.12
    def Poisson(self, rate=1.):

        old_rate = -1

        logfact = np.ones(1024) * -1.
        
        # uses product of uniforms method
        if rate < 5:
            if (rate != old_rate):
                exprate = np.exp(-rate)

            k = -1
            t = 1.
            # trick to doing do-while in Python
            while True:

                k+=1
                t *= self.rand()

                if (t <= exprate):
                    break
        # otherwise uses ratio of uniforms method
        else:
            if (rate != old_rate):
                sqrtrate = np.sqrt(rate)
                lograte = np.log(rate)
            
            # this is a real infinite loop with various break conditions
            while True:
                u = 0.64 * self.rand()
                v = -0.68 + 1.28 * self.rand()

                # outer squeeze for rejection
                if (rate > 13.5):
                    v2 = np.square(v)

                    if (v >=0.):
                        if(v2 > 6.5*u*(0.64-u)*(u+0.2)):
                           continue

                    else:

                        if (v2 > 9.6*u*(0.66-u)*(u+0.07)):
                            continue

                k = int(sqrtrate*(v/u)+rate+0.5)

                if (k < 0):
                        continue

                u2 = np.square(u)
                
                # inner squeeze for acceptance
                if (rate > 13.5):
                    if (v>0.):
                        if (v2 < 15.2*u2*(0.61-u)*(0.8-u)):
                            break

                    else:
                        if (v2 < 6.76*u2*(0.62-u)*(1.4-u)):
                            break
                if (k < 1024):
                    if (logfact[k] < 0.):
                        logfact[k] = special.gammaln(k + 1.)

                    lfac = logfact[k]

                else:
                    lfac = special.gammaln(k+1.)

                # apparently slow, so only does this if absolutely necessary
                p = sqrtrate*np.exp(-rate + k*lograte - lfac)

                if (u2 < p):
                    break

        old_rate = rate
        return k
