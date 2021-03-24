# -*- coding: utf-8 -*-
"""
Class for the results of an optimisation with LAYLA
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

#import sys
#sys.path.append(r'C:\BELLA_and_LAYLA')
#from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

class LAYLA_Results():
    " An object for storing the results of an optimisation with LAYLA"
    def __init__(self, parameters, targets):
        "Initialise the results of an optimisation with LAYLA"
        self.completed = False
        # solution stacking sequence
        self.ss_best = None
        # solution lamination parameters
        self.lampam_best = None
        # solution objective
        self.objective = None
        # stacking sequence solutions at each outer step
        self.ss_tab = np.zeros((
            parameters.n_outer_step, targets.n_plies), int)
        # solution lamination parameters at each outer step
        self.lampam_tab_tab = np.zeros((
            parameters.n_outer_step, 12), float)
        # solution objectives at each outer step
        self.obj_tab = np.NaN*np.ones((
            parameters.n_outer_step,), float)
#        # number of objective function evaluations at each outer step
#        self.n_obj_func_calls_tab = np.NaN*np.ones((
#            parameters.n_outer_step,), int)
        # numbers of stacks at the last level of the last group search
        self.n_designs_last_level_tab = np.NaN*np.ones((
            parameters.n_outer_step,), int)
        # numbers of repaired stacks at the last level of the last group search
        self.n_designs_repaired_tab = np.NaN*np.ones((
            parameters.n_outer_step,), int)
        # numbers of unique repaired stacks at the last group search
        self.n_designs_repaired_unique_tab = np.NaN*np.ones((
            parameters.n_outer_step,), int)
        # numbers of outer steps performed
        self.number_of_outer_steps_performed = None
        # number of the outer step that finds the best solution
        self.n_outer_step_best_solution = None

    def __repr__(self):
        " Display object "

        return f'''
Results with LAYLA:

    Stacking sequence: {self.ss_best}
    Lamination parameters 1-4: {self.lampam_best[:4]}
    Lamination parameters 5-8: {self.lampam_best[4:8]}
    Lamination parameters 9-12: {self.lampam_best[8:]}
'''
