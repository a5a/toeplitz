# -*- coding: utf-8 -*-
"""
DOCTEXT
"""

import numpy


class ToeplitzGP(object):
    # TODO: Use GPy.model as parent here?

    def __init__(self, X, Y, kern):
        # initialise
        raise NotImplementedError

    def set_XY(self, X, Y):
        raise NotImplementedError

    def log_likelihood(self):
        raise NotImplementedError

    def gradient(self):
        raise NotImplementedError

    def train(self):
        raise NotImplementedError

    def predict(self, X_star):
        raise NotImplementedError


