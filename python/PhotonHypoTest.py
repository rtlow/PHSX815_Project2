#! /usr/bin/env python

# imports of external packages to use in our code
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.special as special

#setting matplotlib ticksize
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

#set matplotlib global font size
matplotlib.rcParams['font.size']=14

# import our Random class from python/Random.py file
sys.path.append(".")
from python.MySort import MySort


# helper functions


# given two sets of data and a significance level, plots the histograms
# with significance
def PlotHypotheses(array0, array1, title, alpha, filename='Hypotheses'):

    N0 = len(array0)
    N1 = len(array1)

    hmin = min(array0[0], array1[0])
    hmax = max(array0[N0-1], array1[N1-1])
    
    # normalization weights
    w1 = np.ones_like(array0)/N0
    w2 = np.ones_like(array1)/N1

    fig = plt.figure(figsize=[12,7])
    ax = plt.axes()

    # just 100 bins for visualization
    plt.hist(array0, 100, weights=w1, color='b', alpha=0.5, label='$P(\lambda | H0)$')
    plt.hist(array1, 100, weights=w2, color='g', alpha=0.5, label='$P(\lambda | H1)$')
    
    # since the list is sorted
    # lambda_crit is achieved at 1-alpha percent of the way through array0
    # or at its end
    lambda_crit = array0[min(int((1-alpha)*N0), N0-1)]
    
    # gets the index of the first value in array1 that is above lambda_crit
    first_leftover = np.where( array1 > lambda_crit )[0][0]

    # since the list is sorted
    # beta percent of the way through array1 is the first value
    # in array1 that is above lambda_crit
    # knowing that index, divide by N1 to get beta
    beta = first_leftover/N1

    # fancy plot formatting
    plt.axvline(lambda_crit, color='k')
    plt.text(lambda_crit, ax.get_ylim()[1] * 0.8, '$\\alpha = {:.3f}$'.format(alpha))
    plt.plot([],[], '', label='$\\alpha = {:.3f}$'.format(alpha))
    plt.plot([],[], '', label='$\\beta = {:.3f}$'.format(beta))
    plt.plot([],[], '', label='$\lambda_{crit} = $' + '${:.3f}$'.format(lambda_crit))
    ax.set_yscale('log')

    ax.set_xlabel('$\lambda = \log [ \mathcal{L}(H1) / \mathcal{L}(H0) ]$')
    ax.set_ylabel('Probability')

    plt.legend()

    plt.title(title)
    
    plt.grid(True)

    fig.savefig(filename + '.png')
    plt.show()

    return alpha, beta, lambda_crit

# main function for our CookieAnalysis Python code
if __name__ == "__main__":
   
    haveInput = [False, False]

    InputFile = [None, None]

    alpha = 0.
    
    # reading in the cmd args
    for i in range(1,len(sys.argv)):
        if sys.argv[i] == '-h' or sys.argv[i] == '--help':
            continue

        # seeing if we have input files
        if sys.argv[i] == '-H0':
            InputFile[0] = sys.argv[i+1]
            haveInput[0] = True

        if sys.argv[i] == '-H1':
            InputFile[1] = sys.argv[i+1]
            haveInput[1] = True

        if sys.argv[i] == '-alpha':
            alpha = float(sys.argv[i + 1])
    
    if '-h' in sys.argv or '--help' in sys.argv or not np.all(haveInput):
        print ("Usage: %s [options] -H0 [input file for H0] -H1 [input file for H1]" % sys.argv[0])
        print ("  options:")
        print ("   --help(-h)          print options")
        print ("   -alpha [number]     significance of test")
        print
        sys.exit(1)
    
    # reading in data from files

    counts = [None, None]

    counts[0] = np.loadtxt(InputFile[0])
    counts[1] = np.loadtxt(InputFile[1])

    LLR = []
    Nmeas = 0
    
    # use simulated data for
    # building probability histograms
    
    hists = [None, None]
    bases = [None, None]
    for h in range(2):
        
        # make this into a 1D array
        reshaped = np.reshape( counts[h], -1 )
        
        # making unit-width bins
        bins = np.arange( np.floor( reshaped.min() ), np.ceil( reshaped.max() ) )
        
        # with unit-width bins, the values are true probabilities
        values, base = np.histogram( reshaped, bins=bins, density=True )

        hists[h] = values
        bases[h] = base

    # loop over all hypotheses
    for h in range(2):
        
        this_hyp = []
        
        Nexp = len(counts[h])

        # loop over all experiments
        for e in range(Nexp):
            Nmeas = len(counts[h][e])

            LogLikeRatio = 0.
            
            ok_LLR = True

            # loop over all measurements to calculate the LLR
            for m in range(Nmeas):
                
                prob_of_H0 = 0
                prob_of_H1 = 0

                try:

                    prob_of_H0 = hists[0][np.digitize( counts[h][e][m], bases[0], right=True )]
                    prob_of_H1 = hists[1][np.digitize( counts[h][e][m], bases[1], right=True )]

                    if (prob_of_H0 > 0) and (prob_of_H1 > 0):

                        # LLR is a sum; one contributes positive, other negative
                        LogLikeRatio += np.log( prob_of_H1 ) # LLR for H1

                        LogLikeRatio -= np.log( prob_of_H0 ) # LLR for H0
                    
                    else:
                        #ok_LLR = False
                        continue

                except:

                    #ok_LLR = False
                    continue

            if ok_LLR:

                this_hyp.append(LogLikeRatio)

        LLR.append(this_hyp)
    
    # sort the data
    Sorter = MySort()

    LLR[0] = np.array(Sorter.DefaultSort(LLR[0]))
    LLR[1] = np.array(Sorter.DefaultSort(LLR[1]))

    plot_title = "{} measurements / experiment".format(Nmeas)

    fname = 'signalNoise'
    # plot the histogram
    alpha, beta, lambda_crit = PlotHypotheses(LLR[0], LLR[1], plot_title, alpha, filename=fname)
    
    print
    print('Plot output to file as ' + fname + '.png')
    print('Using alpha = {:.3f}, we get lambda_crit = {:.3f} and beta = {:.3f}'.format(alpha, lambda_crit, beta))
    print
    sys.exit(1)

