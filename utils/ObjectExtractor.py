# -*- coding: utf-8 -*-
"""ObjectExtractor.py

Implements class ObjectExtractor to extract the physics objects and 
their kinematic observables from an uproot tree, opened from some file.
"""

import numpy as np
import pandas as pd

class ObjectExtractor:
    """Class ObjectExtractor

    Extracts all relevant physics properties from
    pandas dataframes, created from an uproot-loaded tree.
    """

    def configure_dfs(self):
        ### Objects with same dimension in dataframe are imported together
        ### I.e. MET, cuts only have 1 entry, muons can have up to 2, jets more
        self.dim1entries_bkg = ['recoPFMetPt', 'recoPFMetPhi', 'cutsVec*', 'genwgt']
        self.dim1entries_sig = ['recoPFMetPt', 'recoPFMetPhi', 'cutsVec*']
        self.dimVertex =  ['recoDr', 'recoVxy', 'recoVz']
        self.dimMu = ['recoPt', 'recoEta', 'recoPhi']
        self.dimJet = ['recoPFJetPt', 'recoPFJetEta', 'recoPFJetPhi']
    
    def __init__(self, uproot_tree, sample_name='', executor=None):
        """Parameters:

            uproot_tree (uproot object): tree loaded from a ROOT file
            sample_name (str): bkg or signal sample name
            executor (Executor): for concurrency (i.e. futures module)
        """
        self.configure_dfs()
        self.sample_name = sample_name
        
        self.uproot_tree = uproot_tree
        self.vertex_df = self.uproot_tree.pandas.df(self.dimVertex, executor=executor)
        self.muons_df = self.uproot_tree.pandas.df(self.dimMu, executor=executor)
        self.jets_df = self.uproot_tree.pandas.df(self.dimJet, executor=executor)
        try:
            self.one_dim_entries_df = self.uproot_tree.pandas.df(self.dim1entries_bkg, executor=executor)
        except KeyError:
            self.one_dim_entries_df = self.uproot_tree.pandas.df(self.dim1entries_sig, executor=executor)
        
    def get_muons(self):
        return self.muons_df.reset_index(level=1)
    
    def get_MET(self):
        return self.one_dim_entries_df[['recoPFMetPt', 'recoPFMetPhi']]
    
    def get_leading_jet(self):
        return self.jets_df.loc[(slice(None),0),slice(None)].reset_index(level=1)
    
    def get_cuts(self):
        return self.one_dim_entries_df.filter(regex='cutsVec')#[[f'cutsVec[{cut}]' for cut in numCuts]]
    
    def get_cuts_crit(self):
        cuts_df = self.get_cuts()
        return [ cuts_df[column] == 1 for column in cuts_df ]
    
    def get_vertex(self):
        return self.vertex_df.reset_index(level=1)
    
    def get_weights(self):
        try:
            pileup_df = self.one_dim_entries_df[['genwgt']] # will also include pileup when available
            genwgt_df = pileup_df['genwgt']
        except KeyError:
            print(f'Sample "{self.sample_name}" does not have either pileup or weight information')
            genwgt_df = pd.Series(np.ones((len(self.one_dim_entries_df))))
            genwgt_df = genwgt_df.rename('genwgt')
            genwgt_df.index.name = 'entry'
        return genwgt_df
    
    def get_all(self):
        objects = {}
        objects['MET'] = self.get_MET()
        objects['cuts'] = self.get_cuts()
        objects['cutsCrit'] = self.get_cuts_crit()
        objects['vertex'] = self.get_vertex()
        objects['muons'] = self.get_muons()
        objects['leadingJet'] = self.get_leading_jet()
        objects['genwgt'] = self.get_weights()
        return objects
