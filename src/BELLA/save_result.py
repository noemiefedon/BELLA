# -*- coding: utf-8 -*-
"""
Function to save the results of multi-panel optimisations

- save_result_BELLAs
    saves the results for the design of a multipanel structure

- save_result_BELLA_one_pdl
    saves the results at one iteration of the design of a multipanel
    structure
"""
import sys
import numpy as np
import pandas as pd

sys.path.append(r'C:\BELLA')
from src.divers.excel import append_df_to_excel
from src.buckling.buckling import buckling_factor
from src.guidelines.ipo_oopo import calc_penalty_ipo_oopo_mp
from src.guidelines.contiguity import calc_penalty_contig_mp
from src.guidelines.disorientation import calc_number_violations_diso_mp
from src.guidelines.ten_percent_rule import calc_penalty_10_ss
from src.CLA.lampam_functions import calc_lampam
from src.CLA.ABD import D_from_lampam

def save_result_BELLAs(filename, multipanel, constraints, parameters,
                       obj_func_param, pdl, output, mat=None, only_best=False):
    """
    saves the results for the design of a multipanel structure
    """
    if only_best:
        table_res = pd.DataFrame()
        if hasattr(output, 'time'):
            table_res.loc[0, 'time (s)'] = output.time
        table_res = save_result_BELLA_one_pdl(
            table_res,
            multipanel,
            constraints,
            parameters,
            obj_func_param,
            output,
            None,
            0,
            mat)
        append_df_to_excel(
            filename, table_res, 'Best Result', index=False, header=True)
        return 0


    for ind in range(parameters.n_ini_ply_drops):
        table_res = pd.DataFrame()
        table_res = save_result_BELLA_one_pdl(
            table_res,
            multipanel,
            constraints,
            parameters,
            obj_func_param,
            output,
            pdl[ind],
            ind,
            mat)
        append_df_to_excel(
            filename, table_res, 'Results', index=False, header=True)
    # to save best result
    table_res = pd.DataFrame()
    table_res.loc[0, 'time (s)'] = output.time
    table_res = save_result_BELLA_one_pdl(
        table_res,
        multipanel,
        constraints,
        parameters,
        obj_func_param,
        output,
        pdl[output.ind_mini],
        0,
        mat)
    append_df_to_excel(
        filename, table_res, 'Best Result', index=False, header=True)
    return 0

def save_result_BELLA_one_pdl(
        table_res, multipanel, constraints, parameters, obj_func_param,
        output, pdl=None, ind=0, mat=None):
    """
    saves the results at one iteration of the design of a multipanel
    structure
    """

    if parameters is None:
        if multipanel.panels[0].N_x == 0 and multipanel.panels[0].N_y == 0:
            save_buckling = False
        else:
            save_buckling = True
    else:
        save_buckling = parameters.save_buckling

    if hasattr(output, 'obj_constraints'):
        table_res.loc[ind, 'obj_constraints'] \
        = output.obj_constraints_tab[ind]
    if hasattr(output, 'n_obj_func_calls_tab'):
        table_res.loc[ind, 'n_obj_func_calls']  \
        = output.n_obj_func_calls_tab[ind]
    # table_res.loc[ind, 'n_designs_last_level'] \
    # = output.n_designs_last_level_tab[ind]
    # table_res.loc[ind, 'n_designs_after_ss_ref_repair'] \
    # = output.n_designs_after_ss_ref_repair_tab[ind]
    # table_res.loc[ind, 'n_designs_after_thick_to_thin'] \
    # = output.n_designs_after_thick_to_thin_tab[ind]
    # table_res.loc[ind, 'n_designs_after_thin_to_thick'] \
    # = output.n_designs_after_thin_to_thick_tab[ind]
    # table_res.loc[ind, 'n_designs_repaired_unique'] \
    # = output.n_designs_repaired_unique_tab[ind]

    table_res.loc[ind, 'penalty_spacing'] = output.penalty_spacing_tab[ind]

    ss = output.ss

    norm_diso_contig = np.array(
        [panel.n_plies for panel in multipanel.panels])

    n_diso = calc_number_violations_diso_mp(ss, constraints)
    penalty_diso = np.zeros((multipanel.n_panels,))
    if constraints.diso and n_diso.any():
        penalty_diso = n_diso / norm_diso_contig
    else:
        penalty_diso = np.zeros((multipanel.n_panels,))

    n_contig = calc_penalty_contig_mp(ss, constraints)
    penalty_contig = np.zeros((multipanel.n_panels,))
    if constraints.contig and n_contig.any():
        penalty_contig = n_contig / norm_diso_contig
    else:
        penalty_contig = np.zeros((multipanel.n_panels,))

    lampam = np.array([calc_lampam(ss[ind_panel]) \
                       for ind_panel in range(multipanel.n_panels)])
        
    penalty_10 = np.zeros((multipanel.n_panels,))
    if constraints.rule_10_percent and constraints.rule_10_Abdalla:
        penalty_10 = calc_penalty_10_ss(ss, constraints, lampam, mp=True)
    else:
        penalty_10 = calc_penalty_10_ss(ss, constraints, LPs=None, mp=True)

    penalty_bal_ipo, penalty_oopo = calc_penalty_ipo_oopo_mp(
        lampam, constraints)

    for ind_p, panel in enumerate(multipanel.panels):

        table_res.loc[ind + ind_p, 'index panel'] = ind_p + 1
        table_res.loc[ind + ind_p, 'n_plies'] = panel.n_plies
        if hasattr(output, 'obj_no_constraints'):
            table_res.loc[ind + ind_p, 'obj_no_constraints'] \
            = output.obj_no_constraints_tab[ind][ind_p]
        table_res.loc[ind + ind_p, 'n_violations_diso'] = n_diso[ind_p]
        table_res.loc[ind + ind_p, 'n_violations_contig'] = n_contig[ind_p]
        table_res.loc[ind + ind_p, 'ipo: |lampam[3]| + |lampam[4]|'] \
        = abs(lampam[ind_p, 2]) + abs(lampam[ind_p, 3])
        table_res.loc[ind + ind_p, 'oopo: |lampam[11]| + |lampam[12]|'] \
        = abs(lampam[ind_p, 10]) + abs(lampam[ind_p, 11])
        table_res.loc[ind + ind_p, 'percentage_0_plies'] \
        = sum(output.sst[ind_p] == 0) / (output.sst[ind_p].size)
        table_res.loc[ind + ind_p, 'percentage_90_plies'] \
        = sum(output.sst[ind_p] == 90) / (output.sst[ind_p].size)
        table_res.loc[ind + ind_p, 'percentage_+45_plies'] \
        = sum(output.sst[ind_p] == 45) / (output.sst[ind_p].size)
        table_res.loc[ind + ind_p, 'percentage_-45_plies'] \
        = sum(output.sst[ind_p] == -45) / (output.sst[ind_p].size)
        table_res.loc[ind + ind_p, 'percentage_+-45_plies'] \
        = (sum(output.sst[ind_p] == 45) + sum(output.sst[ind_p] == -45)) \
        / (output.sst[ind_p].size)

        table_res.loc[ind + ind_p, 'penalty_diso'] = penalty_diso[ind_p]
        table_res.loc[ind + ind_p, 'penalty_contig'] = penalty_contig[ind_p]
        table_res.loc[ind + ind_p, 'penalty_10'] = penalty_10[ind_p]
        table_res.loc[ind + ind_p, 'penalty_ipo'] = penalty_bal_ipo[ind_p]
        table_res.loc[ind + ind_p, 'penalty_oopo'] = penalty_oopo[ind_p]

        if save_buckling:
            table_res.loc[ind + ind_p, 'n_plies_'] = panel.n_plies
            table_res.loc[ind + ind_p, 'lambda buckling'] = buckling_factor(
                lampam=lampam[ind_p],
                mat=mat,
                n_plies=output.ss[ind_p].size,
                N_x=panel.N_x,
                N_y=panel.N_y,
                length_x=panel.length_x,
                length_y=panel.length_y,
                n_modes=10)

        for angle in constraints.set_of_angles:
            table_res.loc[ind + ind_p, 'n_' + str(angle) + 'plies'] \
            = sum(output.sst[ind_p] == angle)

        for ind_lp in range(12):
            table_res.loc[
                ind + ind_p, 'lampam_error[' + str(ind_lp + 1) + ']'] = \
            abs(panel.lampam_target[ind_lp] - lampam[ind_p, ind_lp])

        for ind_lp in range(12):
            table_res.loc[
                ind + ind_p, 'lampam[' + str(ind_lp + 1) + ']'] = \
            lampam[ind_p, ind_lp]
            
        D = D_from_lampam(lampam[ind_p], mat) 
        D_11 = D[0, 0]
        D_22 = D[1, 1]
        D_16 = D[0, 2]
        D_26 = D[1, 2]
        
        table_res.loc[ind + ind_p, 'gamma'] = abs(D_16/(( (D_11**3) * D_22 )**(1/4)))
        table_res.loc[ind + ind_p, 'delta'] = abs(D_26/(( (D_22**3) * D_11 )**(1/4)))
            
        stack = ' '.join(np.array(output.sst[ind_p], dtype=str))
        table_res.loc[
            ind + ind_p, 'stacking sequences with ply drops included'] = stack

        stack = ' '.join(np.array(output.ss[ind_p], dtype=str))
        table_res.loc[ind + ind_p, 'stacking sequences'] = stack

        if pdl is not None:
            stack = ' '.join(np.array(pdl[ind_p]).astype(str))
            table_res.loc[ ind + ind_p, 'initial ply drop layout'] = stack
            

    return table_res
