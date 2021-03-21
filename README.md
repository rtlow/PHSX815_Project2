# PHSX815_Project2

Photon detection is expected to be a Poisson process. However, the rate parameter that noise follows, as well as the rate at which photons fall onto the detector, depend nontrivially on other parameters. The code in this repository simulates the detection of photons on a single pixel for a user-defined number of measurements over a user-defined number of experiments, taking these confounding factors in account.

## Dependencies

- matplotlib
- numpy
- scipy

## Example Usage

The executable files are `python/PhotonCounter.py` and `python/PhotonHypoTest.py`, and can be called from the command line with the -h flag, which will print the options.

The file `python/PhotonCounter.py` generates the Poisson-distributed data according to each model. The user can specify the model using either the `--model0` or `--model1` flag, and can alter the number of Experiments and Measurements per Experiment using the -Nexp and -Nmeas options respectively.

`model0` is the noise model. Noise electrons have an energy spectrum according to the Fermi-Dirac distriubtion. The number of noise electrons detected is the number of electrons above the Fermi level. This number is obtained via numerical integration of the Fermi-Dirac distribution. The `-T` flag allows the user to specify the temperature of Fermi-Dirac distribution.

`model1` is the signal model. Atmospheric seeing causes the observed intensity of a point source to dance around the detector. We model this by sampling the Airy disk using Metropolis-Hastings sampling. The `-Nburn` flag allows the user to specify how many initial samples to exclude. The `-Nskip` flag allows the user to skip a certain number of samples, since the Metropolis-Hastings method may have unwanted periodicity.

The file `python/PhotonHypoTest.py` analyzes two data files from `PhotonCounter.py`. These two data files, with different models, represent two different hypotheses. The user can specify a significance level for hypothesis testing using the `-alpha` option.
