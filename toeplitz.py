# -*- coding: utf-8 -*-
"""
DOCTEXT
"""

import numpy as np


def toeplitz_inverse(m):
    """
    Computes a fast inverse of a Toeplitz matrix (Alg. 1.1) and
    log det m

    mat = numpy.ndarray with Toeplitz form
    """
    c = m[:, 0]
    N = len(c)
    c_inv = np.zeros((N, N))  # intialising the matrix
    v, l = modified_trench_helper(c)

    c_inv[0, :] = v[::-1]
    c_inv[:, 0] = v[::-1]
    c_inv[-1, :] = v
    c_inv[:, -1] = v

    for ii in xrange(1, int(np.floor((N-1)/2)+1)):
        for jj in xrange(ii, N - ii):
            c_inv[ii, jj] = (c_inv[ii-1, jj-1]
                             + (v[N-jj-1] * v[N-ii-1] - v[ii-1] * v[jj-1])
                             / v[-1])
            c_inv[jj, ii] = c_inv[ii, jj]
            c_inv[N-ii-1, N-jj-1] = c_inv[ii, jj]
            c_inv[N-jj-1, N-ii-1] = c_inv[ii, jj]

    return c_inv, l


def modified_trench_helper(c):
    """
    The modified Trench algorithm (Alg. 1.2)

    c is a vector containing the first column of the matrix
    """
    N = len(c)
    wiggle = c[1:]/c[0]
    z, l = modified_durbin(len(wiggle), wiggle)
    l = l + N * np.log(c[0])

    v = np.zeros(N)  # initialising
    v[-1] = 1/((1 + wiggle.dot(z)) * c[0])
    v[:-1] = v[-1] * z[::-1]

    return v, l


def modified_durbin(m, wiggle):
    """
    Modified Durbin algorithm
    (Alg. 1.3)
    """
    z = np.zeros(m)

    z[0] = -wiggle[0]
    beta = 1
    alpha = -wiggle[0]
    l = 0

    for ii in xrange(m-1):
        beta = (1 - alpha**2) * beta
        l = l + np.log(beta)
        if ii == 0:  # there has to be a better way than this...?
            alpha = - (wiggle[ii+1] + wiggle[0] * z[0]) / beta
            z[0] = z[0] + alpha*z[0]
        else:
            # print alpha
            alpha = - (wiggle[ii+1] + wiggle[ii::-1].dot(z[:ii+1])) / beta
            z[:ii+1] = z[:ii+1] + alpha*z[ii::-1]
        z[ii+1] = alpha

    beta = (1 - alpha**2) * beta
    l = l + np.log(beta)

    return z, l


def vec_diag_sum(c):
    """
    Alg. 3.1
    """
    pass


def vec_solve_trace(c):
    """
    Alg. 3.2
    """
    pass
