# -*- coding: utf-8 -*-
"""HistogramContainer.py

Implements class HistogramContainer to hold histogram results computed by
class HistogramCalculator and progressively add statistics.
"""

import numpy as np
import pandas as pd

class HistogramContainer:
    """Class HistogramContainer

    This class keeps track of counts, edges, and weight squared
    for a given observable histogram. It has an add method that
    allows 2 such objects to easily sum together (useful when
    iterating over many files).
    """

    def __init__(self, bins=None, numCuts=np.arange(0,16), weight=1):
        """Parameters:

            bins (int): number of bins for this histogram
            If bins not given, default is 60 and set
            init to false (i.e. object not initialized yet)
            numCuts (numpy array): array with index of cuts
            weight: takes into account xsec, lumi, and relative gen. weight
        """

        self.numCuts = numCuts
        if bins is None:
            self.init = False
            self.bins = 60
        else:
            self.init = True
            self.bins = bins
        self.counts = {}
        self.wgt_sqrd = {}
        self.edges = np.zeros(self.bins+1)
        for cut in self.numCuts:
            self.counts[cut] = np.zeros(self.bins)
            self.wgt_sqrd[cut] = pd.Series([])
        self.weight = weight
    
    def __add__(self, new_hists):
        # Can only add 2 HistogramContainer objects together,
        # or a list of counts, edges, and wgt_squared to a H. C.
        if not isinstance(new_hists, list) and not type(new_hists).__name__ == 'HistogramContainer':
            raise TypeError(f'Trying to add non-list and non-HistogramContainer object'
                            f'to HistogramContainer! Type: {type(new_hists)}')
            
        if type(new_hists).__name__ == 'HistogramContainer':
            new_obj = HistogramContainer(new_hists.bins)
            for cut in self.numCuts:
                new_obj.counts[cut] = new_hists.counts[cut] + self.counts[cut]
                new_obj.wgt_sqrd[cut] = new_hists.wgt_sqrd[cut].add(self.wgt_sqrd[cut], fill_value=0)
            new_obj.edges = new_hists.edges
        elif isinstance(new_hists, list):
            new_obj = HistogramContainer(len(new_hists[0][0]))
            for cut, (counts, edges, wgt_sqrd) in enumerate(new_hists):
                new_obj.counts[cut] = self.counts[cut] + counts
                new_obj.wgt_sqrd[cut] = self.wgt_sqrd[cut].add(wgt_sqrd, fill_value=0)
                if cut == 0:
                    new_obj.edges = edges
        new_obj.calc_max_min()
        return new_obj

    def set_weight(self, weight):
        self.weight = weight
        for cut in self.numCuts:
            self.counts[cut] *= self.weight
            self.wgt_sqrd[cut] *= self.weight**2
    
    def calc_max_min(self):
        self.max = {cut:max(self.counts[cut]) for cut in self.numCuts}
        self.min = {cut:min(self.counts[cut]) for cut in self.numCuts}
        
    def get_max(self):
        self.calc_max_min()
        return self.max
    
    def get_min(self):
        self.calc_max_min()
        return self.min
