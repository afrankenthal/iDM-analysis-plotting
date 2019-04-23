# -*- coding: utf-8 -*-
"""PlotMaker.py

Implements class PlotMaker to make plots from HistogramContainers.
"""

import numpy as np
from skhep.visual import MplPlotter as skh_plt

import utils.HistogramContainer as HC

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class EmptyHistogramError(Error):
    """Exception raised for empty histograms."""
    def __init__(self, message):
        self.message = message

class PlotMaker:
    """Class PlotMaker

    Given some HistogramContainer objects calculated by
    class HistogramCalculator, this class will generate
    the relevant physics plots that we are interested in.
    Two plotting functions are provided, plot_binned_data
    and plot_binned_data_error (preferred). The error version
    includes the sum of weights squared that was calculated
    together with the histograms, and uses the scikit-hep
    visual libraries for plotting errorbars.
    """

    # Default plot options
    def_kwargs = {'density': True, 'log': True, 'histtype':'step'}
    
    def __init__(self, histos):
        """Parameters:

            histos (dict): a dict of dict of HistogramContainer objects
        """

        self.kwargs = PlotMaker.def_kwargs.copy()
        self.histos = histos
        
    def plot_binned_data(self, axis, bin_edges, data, *args, **kwargs):
        kwargs = kwargs + self.kwargs
        #The dataset values are the bin centres
        x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        #The weights are the y-values of the input binned data
        weights = data
        return axis.hist(x, bins=bin_edges, weights=weights, *args, **kwargs)

    def plot_binned_data_error(self, axis, bin_edges, data, wgt_sqrd, *args, **kwargs):
        bin_width = bin_edges[1] - bin_edges[0]
        errors = np.sqrt(wgt_sqrd)
        if 'density' in kwargs and kwargs['density'] == True:
            errors = errors/np.sum(data)/bin_width
        errors = errors.reindex(np.arange(1, len(bin_edges)), fill_value=0)
        #The dataset values are the bin centres
        x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        #The weights are the y-values of the input binned data
        weights = data
        return skh_plt.hist(x, ax=axis, bins=bin_edges, weights=weights, errorbars=errors, *args, **kwargs)
        
    def make_group_plot(self, axis, plot_var, bkg_grps, cut, *args, **kwargs):
        """Plots groups of backgrounds together
        
        Parameters: 
            axis (Axis): which axis to plot on
            plot_var (str): which physics variable to plot
            bkg_grps (dict): whichc groups of bkg samples to plot together
            cut (int): which cut (0-6 usually) to plot
        """

        new_kwargs = {**self.kwargs, **kwargs}
        grp_histos = {}
        for grp_name, bkgs in bkg_grps.items():
            grp_histos[grp_name] = HC.HistogramContainer()
#             if grp == 'QCD': continue
                #or grp == 'V+Jets': continue
            for bkg in bkgs:
                if bkg not in self.histos[plot_var] or not self.histos[plot_var][bkg]: continue
                grp_histos[grp_name] += self.histos[plot_var][bkg]
            if grp_histos[grp_name].init == False or np.sum(grp_histos[grp_name].counts[cut]) == 0:
                grp_histos.pop(grp_name)
        if len(grp_histos) == 0:
            raise EmptyHistogramError(f'No counts for cut {cut}!')
                
        if new_kwargs['density'] == False:
            max_val = 10*max([histo.get_max()[cut] for histo in grp_histos.values()])
            min_val = min([histo.get_min()[cut] for histo in grp_histos.values()])
            min_val = min(0.1*min_val, 1.0)
        else:
            binwidth = next(iter(grp_histos.values())).edges[1] - next(iter(grp_histos.values())).edges[0]
            max_vals = np.array([histo.get_max()[cut]/np.sum(histo.counts[cut])/binwidth for histo in grp_histos.values()])
            min_vals = np.array([histo.get_min()[cut]/np.sum(histo.counts[cut])/binwidth for histo in grp_histos.values()])
            max_val = 10*max(max_vals[~np.isnan(max_vals)])
            min_val = 0.1*min(min_vals[~np.isnan(min_vals)])

        for grp, histo in grp_histos.items():
            if not any(i > 0 for i in histo.counts[cut]): continue
                
            if grp == '60p0' or grp == '5p25':
                new_kwargs['ls'] = '--'
#                 plot_binned_data(axis, edges[grp], counts[grp], label=grp, *args, **kwargs)
            self.plot_binned_data_error(axis, histo.edges, histo.counts[cut], histo.wgt_sqrd[cut],
                                        label=grp, *args, **new_kwargs)
            axis.set_ylim([min_val, max_val])
        
    def make_bkg_plot(self, axis, plot_var, cut, *args, **kwargs):
        """Plots each background sample separately (e.g. QCD_HTXXtoYY)
        
        Parameters: 
            axis (Axis): which axis to plot on
            plot_var (str): which physics variable to plot
            cut (int): which cut (0-6 usually) to plot
        """
        new_kwargs = {**self.kwargs, **kwargs}
        
        if new_kwargs['density'] == False:
            max_val = 10*max([histo.get_max()[cut] for bkg,histo in self.histos[plot_var].items() if bkg != 'DYJetsToTauTau'])
            min_val = min([histo.get_min()[cut] for histo in self.histos[plot_var].values()])
            min_val = min(0.1*min_val, 1.0)
        else:
            binwidth = next(iter(self.histos[plot_var].values())).edges[1] - next(iter(self.histos[plot_var].values())).edges[0]
            max_vals = np.array([histo.get_max()[cut]/np.sum(histo.counts[cut])/binwidth for histo in self.histos[plot_var].values()])
            min_vals = np.array([histo.get_min()[cut]/np.sum(histo.counts[cut])/binwidth for histo in self.histos[plot_var].values()])
            max_val = 10*max(max_vals[~np.isnan(max_vals)])
            min_val = 0.1*min(min_vals[~np.isnan(min_vals)])
        
        for bkg, histo in self.histos[plot_var].items():
            if bkg not in self.histos[plot_var] or not histo: continue
            if bkg == 'DYJetsToTauTau': continue
            plot_binned_data_error(axis, histo.edges, histo.counts[cut], histo.wgt_sqrd[cut],
                                   label=bkg, *args, **new_kwargs)
        axis.set_ylim([min_val, max_val])
