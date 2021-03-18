#! /usr/bin/env python

# imports of external packages to use in our code
import sys
import numpy as np

# import our Random class from python/Random.py file
sys.path.append(".")
from python.Random import Random

# global variables
kB = 8.617333e-5 # [ev/k]

Ef = 1.12 / 2 # band gap of Si divided by 2 [eV]

Nelectrons = 1000 # no. free electrons per pixel

MeanSeeing = 10 # standard deviation of atmospheric seeing [arcsecond]

I0 = 100 # maximum intensity of Airy disk [W/m^2]

aper_a = 0.5 # aperture radius [m]

aper_R = 1 # distance from aperture to focal plane [m] 

lamb = 4000e-10 # mean observation wavelength [Angstrom]

# default seed
seed = 5555


# read the user-provided seed from the command line (if there)
# logic is up here so that random object can be global
if '-seed' in sys.argv:
    p = sys.argv.index('-seed')
    seed = sys.argv[p+1]

# class instance of our Random class using seed
random = Random(seed)


# returns the probability of finding an excited electron
# at energy E [eV] and temperature T [K]
def FermiDirac(E, T):
    return 1/( 1 + np.exp( (E - Ef)/ (kB * T) ) )

# distribution is bounded by 1 on top
def Flat():
    return 1


# provides samples from the Fermi-Dirac distribution
# using hit-miss method
def sampleFermiDirac(Nsample, T):
    
    i = 0

    samples = []
    
    while( i < Nsample ):
        
        # just need to sample over [0, 1]
        X = random.rand()

        R = FermiDirac(X, T)/Flat()

        ran = random.rand()

        # reject the sample and continue
        if (ran > R):
            continue
        else:
            samples.append(X)
            i+=1

    return samples

# 2D Airy Disk
def AiryDisk(x):

    q = np.sqrt( x[0]**2 + x[1]**2 )

    arg = 2 * np.pi * aper_a * q / (lamb * aper_R)

    return I0 * (2 * special.j1(arg) / arg ) ** 2

def Gaussian(x, mu, sig):
    return 1/( sig * np.sqrt(2 * np.pi) ) np.exp( - ( (x - mu)/( 4 * sig ) )**2 )


# TODO Break that logic into its own function so that we can loop easier TODO


# MCMC sampling for the Airy Disk
def sampleAiry(Nburn=0, Nskip=0):
    
    #initial x and y

    x = [1, 1]

    p_x = [ random.Normal(mu=x[0], sig=MeanSeeing) ,\
                  random.Normal(mu=x[1], sig=MeanSeeing) ]

    acceptance_prob = min( 1,\
                          AiryDisk(p_x) * Gaussian(x[0], p_x[0], MeanSeeing) * Gaussian(x[1], p_x[1], MeanSeeing)/ \
                          ( AiryDisk(x) * Gaussian(p_x[0], x[0], MeanSeeing) * Gaussian(p_x[1], x[1], MeanSeeing) ) )
    R = random.rand()



# TODO Write sampling code (MCMC, etc) TODO

# main function for experiment code
if __name__ == "__main__":
    # if the user includes the flag -h or --help print the options
    if '-h' in sys.argv or '--help' in sys.argv:
        print ("Usage: %s [options]" % sys.argv[0] )
        print('Options:')
        print("-seed [number]       random seed") 
        print('-T [number]          temperature in Kelvin')
        print('-Nmeas [number]      no. measurements per experiment')
        print('-Nexp [number]       no. experiments')
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
    if '-output' in sys.argv:
        p = sys.argv.index('-output')
        OutputFileName = sys.argv[p+1]
        doOutputFile = True

    if '--model0' in sys.argv:
        model0 = True

    if '--model1' in sys.argv:
        model1 = True



    if doOutputFile:
        outfile = open(OutputFileName, 'w')
        outfile.write(str(rate)+" \n")
        for e in range(0,Nexp):
            for t in range(0,Nmeas):
                outfile.write(str(random.Poisson(rate))+" ")
            outfile.write(" \n")
        outfile.close()
    else:
        print(rate)
        for e in range(0,Nexp):
            for t in range(0,Nmeas):
                print(random.Poisson(rate), end=' ')
            print(" ")
   
