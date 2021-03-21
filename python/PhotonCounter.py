#! /usr/bin/env python

# imports of external packages to use in our code
import sys
import numpy as np
import scipy.special as special

# import our Random class from python/Random.py file
sys.path.append(".")
from python.Random import Random

# global variables
kB = 8.617333e-5 # [ev/k]

Ef = 1.12 / 2 # band gap of Si divided by 2 [eV]

Nelectrons = 10000 # no. free electrons per pixel

MeanSeeing = 10 # standard deviation of atmospheric seeing [arcsecond]

I0 = 1000 # maximum intensity of Airy disk [counts/second]


# returns the probability of finding an excited electron
# at energy E [eV] and temperature T [K]
def FermiDirac(E, T):
    return 1/( 1 + np.exp( (E - Ef)/ (kB * T) ) )


# number of noise electrons we measure
# are the number of electrons above the
# Fermi level, so we can do numerical
# integration
def sampleFermiDirac(Nsample, T, n=100):
    
    samples = []
    
    for i in range(Nsample):

        # do the integration here
        
        # lower bound is Ef
        a = Ef

        # upper bound is 1

        b = 1

        V = b - a

        int_samples = []

        xi = np.random.uniform(a, b, size=n)

        int_samples.append( FermiDirac(xi, T) )

        integral = (V/n) * np.sum(int_samples)

        rate = Nelectrons * integral

        samp = np.random.poisson(rate)

        samples.append(samp)

    return samples

# 2D Airy Disk
def AiryDisk(x):

    q = np.sqrt( x[0]**2 + x[1]**2 )

    return I0 * (2 * special.j1(q) / q ) ** 2

# just a Gaussian
def Gaussian(x, mu, sig):
    return 1/( sig * np.sqrt(2 * np.pi) ) * np.exp( - ( (x - mu)/( 4 * sig ) )**2 )

# doing the MCMC sampling
def get_MCMC_sample(x):

    p_x = [ np.random.normal(x[0], MeanSeeing) ,\
            np.random.normal(x[1], MeanSeeing) ]
    
    # since sigma is the same, the ratio between
    # gaussian probabilities is just 1
    acceptance_prob = min( 1, AiryDisk(p_x) / AiryDisk(x) )
    R = np.random.uniform()

    if R <= acceptance_prob:
        return p_x
    else:
        return x


# MCMC sampling for the Airy Disk
def sampleAiry(Nsamp, Nburn=0, Nskip=0):
    
    #initial x and y

    x = [1, 1]

    # burning off Nburn samples
    for n in range(Nburn):
        x = get_MCMC_sample(x)
        
    samples = []

    # generating our samples
    for n in range(Nsamp):
        
        # skipping every Nburn value
        for i in range(Nburn):

            x = get_MCMC_sample(x)


        x = get_MCMC_sample(x)
        
        # use the Airy disk to convert these locations to
        # rate parameters
        rate = AiryDisk(x)
        
        # take a sample from Poisson
        samp = np.random.poisson(rate)

        samples.append(samp)

    return samples

# main function for experiment code
if __name__ == "__main__":
    # if the user includes the flag -h or --help print the options
    if '-h' in sys.argv or '--help' in sys.argv:
        print ("Usage: %s [options]" % sys.argv[0] )
        print('Options:')
        print('-T [number]          temperature in Kelvin')
        print('-Nmeas [number]      no. measurements per experiment')
        print('-Nexp [number]       no. experiments')
        print('-Nskip [number]      no. samples to skip in MCMC')
        print('-Nburn [number]      no. samples to burn in MCMC')
        print('-output [string]     output file name')
        print('--model0             simulate model 0')
        print('--model1             simulate model 1')
        print
        sys.exit(1)



    # default number of exposures (letting light collect for fixed time) - per experiment
    Nmeas = 1

    # default number of experiments
    Nexp = 1

    # default temperature
    T = 300

    # output file defaults
    doOutputFile = False

    # do model0 by default
    model0 = True

    # default Nskip
    Nskip = 0

    # default Nburn
    Nburn = 0


    if '-T' in sys.argv:
        p = sys.argv.index('-T')
        ptemp = float(sys.argv[p+1])
        if ptemp > 0:
            T = ptemp
    if '-Nmeas' in sys.argv:
        p = sys.argv.index('-Nmeas')
        Nt = int(sys.argv[p+1])
        if Nt > 0:
            Nmeas = Nt
    if '-Nexp' in sys.argv:
        p = sys.argv.index('-Nexp')
        Ne = int(sys.argv[p+1])
        if Ne > 0:
            Nexp = Ne
    if '-Nskip' in sys.argv:
        p = sys.argv.index('-Nskip')
        Ne = int(sys.argv[p+1])
        if Ne > 0:
            Nskip = Ne
    if '-Nburn' in sys.argv:
        p = sys.argv.index('-Nburn')
        Ne = int(sys.argv[p+1])
        if Ne > 0:
            Nburn = Ne
    if '-output' in sys.argv:
        p = sys.argv.index('-output')
        OutputFileName = sys.argv[p+1]
        doOutputFile = True

    if '--model0' in sys.argv:
        model0 = True

    if '--model1' in sys.argv:
        model0 = False


    experiments = []
    
    for n in range(Nexp):

        if model0:
            
            measurements = sampleFermiDirac(Nmeas, T)
            experiments.append(measurements)

        else:

            measurements = sampleAiry(Nmeas, Nburn, Nskip)
            experiments.append(measurements)

    if doOutputFile:

        np.savetxt(OutputFileName, experiments)


    else:

        print( experiments )
