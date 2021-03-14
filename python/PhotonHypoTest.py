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

# PoisProb approximates the Poisson probability
# using the Stirling approximation
def PoisProb(x, rate):

    logP = x * np.log(rate) - rate - special.gammaln(x + 1)

    P = np.exp(logP)
    
    return P

# given two sets of data and a significance level, plots the histograms
# with significance
def PlotHypotheses(array0, array1, title, alpha, filename='Hypotheses'):

    N0 = len(array0)
    N1 = len(array1)

    hmin = min(array0[0], array1[0])
    hmax = max(array0[N0-1], array1[N1-1])

    fig = plt.figure(figsize=[12,7])
    ax = plt.axes()

    # just 100 bins for visualization
    plt.hist(array0, 100, density=True, color='b', alpha=0.5, label='$P(\lambda | H0)$')
    plt.hist(array1, 100, density=True, color='g', alpha=0.5, label='$P(\lambda | H1)$')
    
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

    Nmeas = 0
    rate = []
    counts = []
    need_rate = True
    
    # loop over all hypotheses (only 2)
    for h in range(2):
        
        need_rate = True
        this_hyp = []
        
        with open(InputFile[h]) as ifile:
            
            # parse each line
            for line in ifile:
                
                # first line is the rate parameter
                if need_rate:
                    need_rate = False
                    rate.append(float(line))
                    continue
            
                # each line is a different experiment
                lineVals = line.split()
                Nmeas = len(lineVals)
                
                this_exp = []
                
                # need to go through all measurements to convert them from string to float
                for m in range(Nmeas):
                    this_exp.append(float(lineVals[m]))
                this_hyp.append(this_exp)

        counts.append(this_hyp)


    LLR = []
    
    # loop over all hypotheses
    for h in range(2):
        
        this_hyp = []

        Nexp = len(counts[h])

        # loop over all experiments
        for e in range(Nexp):
            Nmeas = len(counts[h][e])

            LogLikeRatio = 0.

            # loop over all measurements to calculate the LLR
            for m in range(Nmeas):
    
                # LLR is a sum; one contributes positive, other negative
                LogLikeRatio += np.log( PoisProb( counts[h][e][m], rate[1] ) ) # LLR for H1

                LogLikeRatio -= np.log( PoisProb( counts[h][e][m], rate[0] ) ) # LLR for H0

            this_hyp.append(LogLikeRatio)

        LLR.append(this_hyp)
    
    # sort the data
    Sorter = MySort()

    LLR[0] = np.array(Sorter.DefaultSort(LLR[0]))
    LLR[1] = np.array(Sorter.DefaultSort(LLR[1]))

    plot_title = "{} measurements / experiment with rates $\lambda_0 = {:.2f}$, $\lambda_1 = {:.2f}$ counts / sec".format(Nmeas, rate[0], rate[1])

    fname = 'rate1_{:.2f}rate2_{:.2f}'.format(rate[0], rate[1])
    # plot the histogram
    alpha, beta, lambda_crit = PlotHypotheses(LLR[0], LLR[1], plot_title, alpha, filename=fname)
    
    print
    print('Plot output to file as ' + fname + '.png')
    print('Using alpha = {:.3f}, we get lambda_crit = {:.3f} and beta = {:.3f}'.format(alpha, lambda_crit, beta))
    print
    sys.exit(1)

