import sys
import numpy as np
import scipy
import scipy.stats as st
from scipy.stats import norm
import pandas as pd

from sklearn.gaussian_process import GaussianProcess
import pdb
from gaussian_periodic import *


#############################################################################
# COMPUTE_EI.py 

# Computes expected 

# Celestial object models for the following types:
# 	- variable stars = Gaussian Process
#	- supernovae = template

# author: Anita Mehrotra, anitamehrotra@fas.harvard.edu
# date: January 21, 2014
# department: CSE SEAS, Harvard University
# project: CHILE 2014
#############################################################################


# covert Lomb-Scargle significance to proportions
def convert_sig(sig):
    """
    Parameters
    ----------
    sig: a numpy array containing float64's
    """
    return sig/sum(sig)

# compute Phat
def Phat(d):
    """
    Parameters
    ----------
    d: a multidimensional array containing values for a particular 
        theta_j at different times e
    """
    return np.mean(d, axis=1)

# create full list of thetas, uncertainties
def full_list(*args):
    """
    Parameters
    ----------
    items: list of items that need to be expanded
    prob: associated percents or probabilities (e.g. [.5,.1,.4])
    size_full_list: the size of the ultimately expanded list (e.g. 100)

    Output
    ----------
    An expanded list of items proportional to the percentages inputted, e.g.
        [item1, item1, ..., item1, item2, item2, ..., item2, item3, ..item3]
        where number of item1's = e.g. 50, number of item2's = e.g. 10, number
        of item3's = e.g. 40
    """
    items, prob, size_full_list = args[0], args[1], args[2]
    full_list = []
    curr_index = 0
    vals = np.floor(prob*size_full_list)
    
    for i in xrange(len(vals)):
        curr_index = curr_index + vals[i]
        intermediate = np.ones(int(vals[i]))*items[i]
        full_list.append(intermediate)
        if i==len(vals)-1:
            intermediate = np.ones(int(size_full_list-curr_index))*items[i]
            full_list.append(intermediate)
    
    return np.hstack(full_list)


def main():

    # load data
    data_input_file = sys.argv[1]
    theta_input_file = sys.argv[2]

    #data_input_file  = '/Users/anita/Documents/CHILE/data/test_lc/Obj133_148.304044_1.378136'
    #theta_input_file = '/Users/anita/Documents/CHILE/data/test_theta.txt'

    jd, mag, mag_err = np.loadtxt(data_input_file, unpack=True, usecols=[0, 1, 2])
    obj_id, theta0, theta1, theta2, dp0, dp1, dp2, sig0, sig1, sig2 = np.loadtxt(theta_input_file, unpack=True, usecols=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    

    ##### this portion for testing only #####
    first = 16
    jd = jd[:first]
    mag = mag[:first]
    mag_err = mag_err[:first]
    ##### end testing portion #####

    # initial vars
    N = 50  # number of thetas i.e. i = {1, 2, ..., N}
    M = 100 # number of new data points d_j i.e. j = {1, 2, ..., M}
    times = np.arange(0, 3, 0.1)
    len_times = len(times)
    ptildes = np.zeros([N, len_times])

    # work with one object from the second file for now
    obj_id = obj_id[0]
    theta0 = theta0[0]
    theta1 = theta1[0]
    theta2 = theta2[0]
    dp0 = dp0[0]
    dp1 = dp1[0]
    dp2 = dp2[0]
    sig0 = sig0[0]
    sig1 = sig1[0]
    sig2 = sig2[0]

    sig = np.array([sig0, sig1, sig2], dtype=float)
    uncert = np.array([dp0, dp1, dp2], dtype=float)
    thetas = np.array([theta0, theta1, theta2], dtype=float)

    # compute probabilities 
    prob = sig/sum(sig)

    # produce full lists; needed to run Gaussian Process
    allthetas = full_list(thetas, prob, N)
    alluncertainty = full_list(uncert, prob, N)

    # run Gaussian Process
    for i in xrange(N):
        print "i =", i # track progress
        period = allthetas[i]
        period_uncertainty = alluncertainty[i]
        
        gp = GaussianPeriodic(jd, mag, mag_err, period, period_uncertainty, verbose=True)
        querytimes = jd[-1] + (jd[-1] - jd[0])*times
        
        # produce "new data" (samples from Gaussian Process)
        dj = gp.sample_d(M,querytimes)

        # compute phat per theta
        ptildes[i,:] = Phat(dj)

    # compute Expected Information per e and output to csv
    print "shape(ptildes)", np.shape(ptildes)
    lnptildes = np.log(ptildes)
    print "shape(lnptildes)", np.shape(lnptildes)
    EI_e = np.mean(lnptildes, axis=0)
    print "len(lnptildes)", len(lnptildes)
    np.savetxt("expected_information.csv", EI_e, delimiter=',')
    print "len_times", len_times
    print "Finished! See expected_information.csv for output."

    return 0


if __name__ == '__main__':
    main()




