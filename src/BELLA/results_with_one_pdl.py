# -*- coding: utf-8 -*-
"""
Class for the results of an optimisation with BELLA
with one initial ply-drop layout
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

#import sys
#sys.path.append(r'C:\BELLA')
#from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

class BELLA_Results():
    " An object for storing the results of an optimisation with BELLA"

    def __init__(self, parameters, constraints, multipanel):
        "Initialise the results of an optimisation with BELLA"
        self.obj_constraints_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops,), dtype=float)
        self.obj_no_constraints_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops, multipanel.n_panels), dtype=float)

        self.penalty_contig_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops, multipanel.n_panels), dtype=float)

        self.penalty_diso_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops, multipanel.n_panels), dtype=float)

        self.penalty_10_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops, multipanel.n_panels), dtype=float)

        self.penalty_bal_ipo_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops, multipanel.n_panels), dtype=float)

        self.penalty_oopo_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops, multipanel.n_panels), dtype=float)

        self.n_contig_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops, multipanel.n_panels), dtype=int)

        self.n_diso_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops, multipanel.n_panels), dtype=int)

        self.n_obj_func_calls_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops,), int)

        self.n_designs_last_level_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops,), int)

        self.n_designs_after_ss_ref_repair_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops,), int)

        self.n_designs_after_thick_to_thin_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops,), int)

        self.n_designs_after_thin_to_thick_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops,), int)

        self.n_designs_repaired_unique_tab = np.NaN*np.ones((
            parameters.n_ini_ply_drops,), int)

        self.lampam_tab_tab = np.zeros((
            multipanel.n_panels, parameters.n_ini_ply_drops, 12), float)

        self.n_plies_per_angle_tab = np.zeros((
            parameters.n_ini_ply_drops, multipanel.n_panels,
            constraints.n_set_of_angles), float)

        # Initialisation of the array storing all the best stacking sequence
        # solutions: ss_void
        ss_void = []
        for panel in multipanel.panels:
            ss_void.append(np.zeros((panel.n_plies,), dtype=int))
        # Initialisation of the array storing all the stacking sequence solutions:
        # ss_tab
        self.ss_tab = [[]]*(parameters.n_ini_ply_drops)
        for outer_step in range(parameters.n_ini_ply_drops):
            self.ss_tab[outer_step] = ss_void
        # Initialisation of the array storing all the stacking sequence tables:
        # ss_tab_tab
        if constraints.sym \
        and multipanel.n_plies_max % 2 == 0 \
        and sum([p.middle_ply for p in multipanel.panels]) != 0:
            self.ss_tab_tab = np.zeros((
                parameters.n_ini_ply_drops,
                multipanel.n_panels,
                multipanel.n_plies_max + 1), dtype=int)
        else:
            self.ss_tab_tab = np.zeros((
                parameters.n_ini_ply_drops,
                multipanel.n_panels,
                multipanel.n_plies_max), dtype=int)

    def update(self, outer_step, results_one_pdl):
        "Update the results from an optimisation with one ply-drop layout"
        self.ss_tab[outer_step] = results_one_pdl.ss
        self.ss_tab_tab[outer_step] = results_one_pdl.sst

        self.lampam_tab_tab[:, outer_step, :] = results_one_pdl.lampam

        self.obj_constraints_tab[
            outer_step] = results_one_pdl.obj_constraints
        self.obj_no_constraints_tab[
            outer_step] = results_one_pdl.obj_no_constraints

        self.penalty_diso_tab[
            outer_step] = results_one_pdl.penalty_diso
        self.penalty_contig_tab[
            outer_step] = results_one_pdl.penalty_contig
        self.penalty_10_tab[
            outer_step] = results_one_pdl.penalty_10
        self.penalty_bal_ipo_tab[
            outer_step] = results_one_pdl.penalty_bal_ipo
        self.penalty_oopo_tab[
            outer_step] = results_one_pdl.penalty_oopo
        self.n_diso_tab[outer_step] = results_one_pdl.n_diso
        self.n_contig_tab[outer_step] = results_one_pdl.n_contig

        self.n_plies_per_angle_tab[
            outer_step] = results_one_pdl.n_plies_per_angle

        self.n_obj_func_calls_tab[
            outer_step] = results_one_pdl.n_obj_func_calls
        self.n_designs_last_level_tab[
            outer_step] = results_one_pdl.n_designs_last_level
        self.n_designs_after_ss_ref_repair_tab[
            outer_step] = results_one_pdl.n_designs_after_ss_ref_repair
        self.n_designs_after_thick_to_thin_tab[
            outer_step] = results_one_pdl.n_designs_after_thick_to_thin
        self.n_designs_after_thin_to_thick_tab[
            outer_step] = results_one_pdl.n_designs_after_thin_to_thick
        self.n_designs_repaired_unique_tab[
            outer_step] = results_one_pdl.n_designs_repaired_unique

    def __repr__(self):
        " Display object "

        return f'''
Results with BELLA:

***
'''
