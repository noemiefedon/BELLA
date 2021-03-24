# -*- coding: utf-8 -*-
"""
Class for optimisation targets
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

class Targets():
    """
    An object for storing the targets of laminate lay-up optimisations
    """
    def __init__(self, n_plies, lampam=np.zeros((12,), float), stack=None):
        " Create targets of laminate lay-up optimisations"
        self.lampam_initial = np.copy(lampam)
        self.lampam = np.copy(lampam)
        self.n_plies = n_plies
        self.stack = stack
        
    def filter_target_lampams(self, constraints):
        """
        filters applied to the lamination parameters to account for orthotropy 
        """
        # If symmetry is desired, the corresponding target amination parameters 
        # must be set to 0
        if constraints.sym:
            self.lampam[4:8] = 0
        # If the in-plane orthotropy is desired, the corresponding target
        # lamination parameters must be set to 0
        if constraints.ipo:
            self.lampam[2] = 0
            self.lampam[3] = 0
        # If the out-of-plane orthotropy is desired, the corresponding target
        # lamination parameters must be set to 0
        if constraints.oopo:
            self.lampam[10] = 0
            self.lampam[11] = 0
        
    def __repr__(self):
        " Display object "
        return f'''
Targets:
    
    Lamination parameters 1-4: {self.lampam[:4]}
    Lamination parameters 5-8: {self.lampam[4:8]}
    Lamination parameters 9-12: {self.lampam[8:]}
    Number of plies: {self.n_plies}
'''
