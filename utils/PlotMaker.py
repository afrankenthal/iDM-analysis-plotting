# -*- coding: utf-8 -*-
"""PlotMaker.py

Implements class PlotMaker to make plots from HistogramContainers.
"""

import numpy as np
from skhep.visual import MplPlotter as skh_plt
from collections import OrderedDict

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
    
    def __init__(self, histos, bkgs, signals=[]):
        """Parameters:

            histos (dict): a dict of dict of HistogramContainer objects
        """

        self.kwargs = PlotMaker.def_kwargs.copy()
        self.histos = histos
        self.bkgs = bkgs
        self.mchis = signals
        
    def plot_binned_data(self, axis, bin_edges, data, *args, **kwargs):
        kwargs = kwargs + self.kwargs
        #The dataset values are the bin centres
        x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        #The weights are the y-values of the input binned data
        weights = data
        return axis.hist(x, bins=bin_edges, weights=weights, *args, **kwargs)

    def plot_binned_data_error(self, axis, bin_edges, data, wgt_sqrd, *args, **kwargs):
        binwidth = bin_edges[1] - bin_edges[0]
        errors = np.sqrt(wgt_sqrd)
        if 'density' in kwargs and kwargs['density'] == True:
            errors = errors/np.sum(data)/binwidth
        errors = errors.reindex(np.arange(1, len(bin_edges)), fill_value=0)
        #The dataset values are the bin centres
        x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        #The weights are the y-values of the input binned data
        weights = data
        return skh_plt.hist(x, ax=axis, bins=bin_edges, weights=weights, errorbars=errors, *args, **kwargs)

    def plot_stacked_binned_data_error(self, axis, bin_edges, data, wgt_sqrd, *args, **kwargs):
        errors = wgt_sqrd[0]
        for i in np.arange(1, len(wgt_sqrd)):
            errors = errors.add(wgt_sqrd[i], fill_value=0)
        errors = np.sqrt(errors)
        errors = np.array(errors.reindex(np.arange(1, len(bin_edges)), fill_value=0))
        #The dataset values are the bin centres
        x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        x = np.array([x]).repeat(len(data), axis=0)
        x = np.transpose(x)
        #The weights are the y-values of the input binned data
        weights = np.transpose(data)
        return skh_plt.hist(x, ax=axis, bins=bin_edges, weights=weights, errorbars=errors, stacked=True, *args, **kwargs)

    def make_group_stacked_plot(self, axis, plot_var, cut, *args, **kwargs):
        grp_histos = {}
        for bkg, properties in self.bkgs.items():
            grp = properties['group']
            if grp not in grp_histos:
                grp_histos[grp] = HC.HistogramContainer()
            # self.histos[plot_var][bkg].set_weight(properties['weight'])
            # FIXME placeholder while H.C. doesn't have set_weight
            grp_histos[grp].counts[cut] += self.histos[plot_var][bkg].counts[cut] * properties['weight']
            grp_histos[grp].edges = self.histos[plot_var][bkg].edges
            grp_histos[grp].wgt_sqrd[cut] = grp_histos[grp].wgt_sqrd[cut].add(self.histos[plot_var][bkg].wgt_sqrd[cut] * properties['weight']**2, fill_value=0)

        labels = []
        sorted_grp_histos = OrderedDict()
        sorted_keys = sorted(grp_histos, key=lambda obj: max(grp_histos[obj].counts[cut]))
        for key in sorted_keys:
            sorted_grp_histos[key] = grp_histos[key]
            labels.append(key)

        # Log doesn't seem to work with stacked option in scikit-hep
        # so we set it ourselves, nbd
        if 'log' in kwargs and kwargs['log'] == True:
            kwargs.pop('log')
            axis.set_yscale('log', nonposy='clip')

        if 'density' in kwargs and kwargs['density'] == True:
            binwidth = next(iter(sorted_grp_histos.values())).edges[1] - next(iter(sorted_grp_histos.values())).edges[0]
            max_vals = np.array([histo.get_max()[cut]/np.sum(histo.counts[cut])/binwidth for histo in sorted_grp_histos.values()])
            min_vals = np.array([histo.get_min()[cut]/np.sum(histo.counts[cut])/binwidth for histo in sorted_grp_histos.values()])
            max_val = 10*max(max_vals[~np.isnan(max_vals)])
            min_val = 0.1*min(min_vals[~np.isnan(min_vals)])
        else:
            max_val = 10*max([histo.get_max()[cut] for histo in sorted_grp_histos.values()])
            min_val = min([histo.get_min()[cut] for histo in sorted_grp_histos.values()])
            min_val = max(min(0.1*min_val, 1.0), 0.1)

        axis.set_ylim([min_val, max_val])

        self.plot_stacked_binned_data_error(
                axis, next(iter(self.histos[plot_var].values())).edges,
                np.array([sorted_grp_histos[grp].counts[cut] for grp in sorted_grp_histos]),
                [sorted_grp_histos[grp].wgt_sqrd[cut] for grp in sorted_grp_histos], label=labels, *args, **kwargs
                )
        
    def make_group_plot(self, axis, plot_var, cut, *args, **kwargs):
        """Plots groups of backgrounds together
        
        Parameters: 
            axis (Axis): axis to plot on
            plot_var (str): which physics variable to plot
            cut (int): which cut (0-5 usually) to plot
        """

        new_kwargs = {**self.kwargs, **kwargs}
        grp_histos = {}
        for bkg, properties in self.bkgs.items():
            grp = properties['group']
            if grp not in grp_histos:
                grp_histos[grp] = HC.HistogramContainer()
            # self.histos[plot_var][bkg].set_weight(properties['weight'])
            # FIXME placeholder while H.C. doesn't have set_weight
            grp_histos[grp].counts[cut] += self.histos[plot_var][bkg].counts[cut] * properties['weight']
            grp_histos[grp].edges = self.histos[plot_var][bkg].edges
            grp_histos[grp].wgt_sqrd[cut] += self.histos[plot_var][bkg].wgt_sqrd[cut] * properties['weight']**2
            # grp_histos[grp] += self.histos[plot_var][bkg]

        for mchi in self.mchis:
            grp_histos[mchi] = HC.HistogramContainer()
            grp_histos[mchi] += self.histos[plot_var][mchi]
                
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
            if grp in self.mchis:
                new_kwargs['ls'] = ':'
            #plot_binned_data(axis, edges[grp], counts[grp], label=grp, *args, **kwargs)
            self.plot_binned_data_error(axis, histo.edges, histo.counts[cut], histo.wgt_sqrd[cut], label=grp, *args, **new_kwargs)
        axis.set_ylim([min_val, max_val])
        
    def make_bkg_plot(self, axis, plot_var, cut, *args, **kwargs):
        """Plots each background sample separately (e.g. QCD_HTXXtoYY)
        
        Parameters: 
            axis (Axis): axis to plot on
            plot_var (str): which physics variable to plot
            cut (int): which cut (0-5 usually) to plot
        """
        new_kwargs = {**self.kwargs, **kwargs}
        
        if new_kwargs['density'] == False:
            max_val = 10*max([histo.get_max()[cut] for bkg,histo in self.histos[plot_var].items()])
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
            plot_binned_data_error(axis, histo.edges, histo.counts[cut], histo.wgt_sqrd[cut], label=bkg, *args, **new_kwargs)
        axis.set_ylim([min_val, max_val])
