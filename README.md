# PHSX815_Project2

Photon detection is expected to be a Poisson process. However, the rate parameter that noise follows, as well as the rate at which photons fall onto the detector, depend nontrivially on other parameters. The code in this repository simulates the detection of photons on a single pixel for a user-defined number of measurements over a user-defined number of experiments, taking these confounding factors in account.

## Dependencies

- matplotlib
- numpy
- scipy

## Example Usage

The executable files are `python/PhotonCounter.py` and `python/PhotonHypoTest.py`, and can be called from the command line with the -h flag, which will print the options.

The file `python/PhotonCounter.py` generates the Poisson-distributed data. The user can provide a Poisson rate parameter using the `-rate` option, and can alter the number of Experiments and Measurements per Experiment using the -Nexp and -Nmeas options respectively.

The file `python/PhotonHypoTest.py` analyzes two data files from `PhotonCounter.py`. These two data files, with different rate parameters, represent two different hypotheses. The user can specify a significance level for hypothesis testing using the `-alpha` option.
