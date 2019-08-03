# -*- coding: utf-8 -*-
"""HistogramCalculator.py 

Implements HistogramCalculator class and some helper functions
to compute more complex observables with pandas apply method.
"""

import math
import numpy as np
import pandas as pd
import operator
from functools import reduce
import multiprocessing
import concurrent.futures

num_cores = multiprocessing.cpu_count()

# Helper functions to calculate average angles
# This takes a few seconds to run, since we 
# are using the apply method

def parallelize(data, func):
#     data_split = np.array_split(data, partitions)
    pool = multiprocessing.Pool(int(num_cores/2))
#     data = pd.concat(pool.map(func, data_split))
    data = pd.concat(pool.map(func, [group for name, group in data]))
    pool.close()
    pool.join()
    return data

def calcAvgAngle(group):
    # FIXME need to ensure at least 2 muons (otherwise index -1 == 0)
    x = np.cos(group['reco_mu_phi'].iloc[0]) + np.cos(group['reco_mu_phi'].iloc[-1])
    y = np.sin(group['reco_mu_phi'].iloc[0]) + np.sin(group['reco_mu_phi'].iloc[-1])
    return math.atan2(y/2, x/2)

def func_group_apply(df):
    # Applies above function on event-by-event basis
    return df.groupby('entry').apply(calcAvgAngle)

def reducephi(row):
    # Helper function to normalize angle differences to [-Pi, Pi]
    # cond: if abs(phidiff) > Pi => phidiff = phidiff - 2*Pi*sign(phidiff)
    if abs(row) > math.pi:
        return row - 2*math.pi*(row/abs(row))
    return row

class HistogramCalculator:
    """Class HistogramCalculator -- computes physics histograms

    Aggregrates the data into histograms, computing the relevant
    observables. For each file in background or signal samples,
    the class HistogramCalculator is called to compute the histograms.
    """

    def __init__(self, objects, sample_name, sample_type='', numCuts=np.arange(0,16)):
        """Parameters:
            objects (dict): dict of pandas dataframes
            sample_name (str): which bkg or signal sample
            sample_type (str): 'bkg' or 'signal'
            numCuts (numpy array): array of cut indexes
        """
        
        self.sample_name = sample_name
        self.sample_type = sample_type
        self.objects = objects
        self.numCuts = numCuts
        self.cuts = objects['cuts']
        self.cuts_crit = objects['cutsCrit']
        self.MET = objects['MET']
        self.muons = objects['muons']
        self.jet = objects['leadingJet']
        self.vertex = objects['vertex'].reset_index()
        try:
            self.genwgt = objects['gen_wgt']
        except KeyError:
            self.genwgt = pd.Series(np.ones(len(self.cuts)))
            self.genwgt = self.genwgt.rename('gen_wgt')
            self.genwgt.index.name = 'entry'
            
    def cutflows(self):
        incl = np.zeros(len(self.numCuts))
        excl = np.zeros(len(self.numCuts))
        for cut in self.numCuts:
            # For inclusive, apply boolean '&&' with all cuts from 0 to cut
            cuts_to_apply = reduce(operator.and_, self.cuts_crit[0:cut+1])
            incl[cut] = len(self.cuts[cuts_to_apply])
            # For exclusive, apply each cut separately
            cuts_to_apply = self.cuts_crit[cut]
            excl[cut] = len(self.cuts[cuts_to_apply])
        return (incl, excl)
        
    def compute_hist(self, variable_df, **kwargs):
        # Given a dataframe for some observable, adds the
        # gen weight and computes the histogram for it
        if 'bins' not in kwargs:
            kwargs['bins'] = 60
        temp_df = pd.concat([variable_df, self.genwgt],axis=1).dropna() 
        #temp_df = pd.concat([variable_df, self.genwgt], axis=1).dropna()
        temp_df['genwgt_sqrd'] = temp_df['gen_wgt']**2
        counts = {}; edges = {}; wgt_sqrd = {}
        for cut in self.numCuts:
            cuts_to_apply = slice(None) if self.cuts_crit is None else reduce(operator.and_, self.cuts_crit[0:cut+1])
            kwargs['weights'] = temp_df[cuts_to_apply]['gen_wgt']
            counts[cut], edges[cut] = np.histogram(temp_df[cuts_to_apply][variable_df.name], **kwargs)
            # Digitizes data to find out which bin of histogram each row falls in
            bin_idxs = np.digitize(temp_df[cuts_to_apply][variable_df.name], edges[cut])
            temp_df['bin_idx'] = pd.Series(bin_idxs)
            # Uses indexes from above to sum the gen weights squared (for errors)
            wgt_sqrd[cut] = np.sum(temp_df.groupby('bin_idx'))['genwgt_sqrd']
        return list(zip(counts.values(), edges.values(), wgt_sqrd.values()))
    
    def metmuphi(self):
        # Divide data into 'chunks' to more efficiently use parallelization
        muons = self.muons.reset_index()
        muons['data_chunk'] = muons['entry'].mod(int(num_cores * 3 / 2)) # num_cores/2 * 3 chunks/core
        muons = muons.set_index(['entry'])
        # Here, group by data_chunk instead of entry, inside func_group_apply 
        # we also have a groupby('entry')
        avg_muon_angle = parallelize(muons.groupby('data_chunk'), func_group_apply)
        angle_diff = (self.MET['reco_PF_MET_phi'].dropna() - avg_muon_angle).dropna()
        reduced_angle_diff = angle_diff.apply(reducephi).dropna()
        reduced_angle_diff.name = 'reducedAngleDiff'
        return self.compute_hist(reduced_angle_diff, range=(-math.pi, math.pi))
    
    def metjetphi(self):
        angle = (self.MET['reco_PF_MET_phi'] - self.jet['reco_PF_jet_phi']).dropna()
        reduced_angle = angle.apply(reducephi)
        reduced_angle.name = 'reducedAngle'
        return self.compute_hist(reduced_angle, range=(-math.pi, math.pi))
    
    def metpt(self):
        return self.compute_hist(self.MET['reco_PF_MET_pt'], range=(0,2500))
    
    def jetpt(self):
        return self.compute_hist(self.jet['reco_PF_jet_pt'], range=(0,2500))
    
    def leadingmupt(self):
        return self.compute_hist(self.muons['reco_mu_pt'].groupby('entry').nth(0), range=(0,700))

    def subleadingmupt(self):
        return self.compute_hist(self.muons['reco_mu_pt'].groupby('entry').nth(1), range=(0,700))

    def recodr(self):
        return self.compute_hist(self.vertex['reco_vertex_dR'], range=(0,2*math.pi))

    def recovertex(self):
        vertex = np.sqrt(self.vertex['reco_vertex_vxy']**2 + self.vertex['reco_vertex_vz']**2)
        vertex.name = 'vertex'
        return self.compute_hist(vertex, range=(0,300))
