# -*- coding: utf-8 -*-
"""
Function to save laminate design set-up

- save_objective_function_BELLA:
    saves the objective function parameters on Sheet [Objective function]

- save_multipanel:
    saves the data of the multipanel structure:
        - panel geometry
        - panel thickness targets
        - panel lamination parameter targets
        - lamination parameter first-level sensitivities
        - boundaries accross panels

- save_constraints_BELLA
    save the design and manufacturing constraints on Sheet [Constraints]

- save_parameters_BELLA
    saves the optimiser parameters on Sheet [Parameters]

- save_materials
    saves the material properties on Sheet [Materials]
"""
import sys
import numpy as np
import pandas as pd

sys.path.append(r'C:\BELLA')
from src.divers.excel import append_df_to_excel
from src.CLA.lampam_functions import calc_lampam
from src.BELLA.format_pdl import convert_sst_to_ss
from src.guidelines.ipo_oopo import calc_penalty_ipo_oopo_mp
from src.guidelines.contiguity import calc_penalty_contig_mp
from src.guidelines.disorientation import calc_number_violations_diso_mp
from src.guidelines.ten_percent_rule import calc_penalty_10_ss
from src.guidelines.ply_drop_spacing import calc_penalty_spacing
from src.buckling.buckling import buckling_factor

def save_materials(filename, materials):
    """
    saves the material properties on Sheet [Materials]
    """
    table_mat = pd.DataFrame()
    table_mat.loc[0, 'E11'] = materials.E11
    table_mat.loc[0, 'E22'] = materials.E22
    table_mat.loc[0, 'G12'] = materials.G12
    table_mat.loc[0, 'nu12'] = materials.nu12
    table_mat.loc[0, 'nu21'] = materials.nu21
    table_mat.loc[0, 'areal density'] = materials.density_area
    table_mat.loc[0, 'volumic density'] = materials.density_volume
    table_mat.loc[0, 'ply thickness'] = materials.ply_t

    table_mat.loc[0, 'Q11'] = materials.Q11
    table_mat.loc[0, 'Q12'] = materials.Q12
    table_mat.loc[0, 'Q22'] = materials.Q22
    table_mat.loc[0, 'Q66'] = materials.Q66

    table_mat.loc[0, 'U1'] = materials.U1
    table_mat.loc[0, 'U2'] = materials.U2
    table_mat.loc[0, 'U3'] = materials.U3
    table_mat.loc[0, 'U4'] = materials.U4
    table_mat.loc[0, 'U5'] = materials.U5

    table_mat = table_mat.transpose()

    append_df_to_excel(
        filename, table_mat, 'Materials', index=True, header=False)


def save_multipanel(
        filename, multipanel, obj_func_param, sst=None,
        calc_penalties=False, constraints=None, mat=None, save_buckling=False):
    """
    saves the data of the multipanel structure:
        - panel geometry
        - panel thickness targets
        - panel lamination-parameter targets
        - lamination parameter first-level sensitivities
        - boundaries accross panels
        - constraints: design guidelines
        - sst: stacking sequence table
    """
    table_mp = pd.DataFrame()
    table_mp.loc[0, 'Number of panels'] = multipanel.n_panels
    table_mp.loc[0, 'Number of plies max'] = multipanel.n_plies_max
    table_mp.loc[0, 'Area'] = multipanel.area_patches
    table_mp.loc[0, 'Area of all patches'] = multipanel.area_patches
    if mat is not None:
        table_mp.loc[0, 'Weight'] = multipanel.calc_weight(mat.density_area)
    table_mp.loc[0, 'Index of one thickest panel'] = multipanel.ind_thick
    table_mp.loc[0, 'Number of plies max'] = multipanel.n_plies_max

    if calc_penalties:
        # penalty_spacing
        penalty_spacing = calc_penalty_spacing(
            pdl=sst,
            multipanel=multipanel,
            constraints=constraints,
            on_blending_strip=False)
        table_mp.loc[0, 'Penalty spacing'] = penalty_spacing

    table_mp = table_mp.transpose()

    append_df_to_excel(
        filename, table_mp, 'Multipanel', index=True, header=False)

    table_p = pd.DataFrame()

    for ind_p, panel in enumerate(multipanel.panels):

        table_p.loc[ind_p, 'Panel ID'] = panel.ID
        table_p.loc[ind_p, 'Neighbour panel IDs'] \
        = " ".join(np.array(panel.neighbour_panels).astype(str))
        table_p.loc[ind_p, 'Number of plies'] = panel.n_plies

        table_p.loc[ind_p, 'Weighting in MP objective funtion'] = panel.weighting

        if panel.length_x and panel.length_y:
            table_p.loc[ind_p, 'Length_x'] = panel.length_x
            table_p.loc[ind_p, 'Length_y'] = panel.length_y
            table_p.loc[ind_p, 'Area'] = panel.area
        else:
            table_p.loc[ind_p, 'Area'] = panel.area

        if hasattr(panel, 'N_x'):
            table_p.loc[ind_p, 'N_x'] = panel.N_x
            table_p.loc[ind_p, 'N_y'] = panel.N_y

        if hasattr(panel, 'Weight'):
            table_p.loc[ind_p, 'Weight'] = panel.calc_weight(mat.density_area)

        for ind in range(12):
            table_p.loc[ind_p, 'lampam_target[' + str(ind + 1) + ']'] \
            = panel.lampam_target[ind]
        for ind in range(12):
            table_p.loc[ind_p, 'lampam_weightings_ini[' + str(ind + 1) + ']'] \
            = panel.lampam_weightings_ini[ind]
        for ind in range(12):
            table_p.loc[ind_p, 'lampam_weightings[' + str(ind + 1) + ']'] \
            = panel.lampam_weightings[ind]

    if calc_penalties:
        ss = np.array(convert_sst_to_ss(sst))

        norm_diso_contig = np.array(
                [panel.n_plies for panel in multipanel.panels])
        n_diso = calc_number_violations_diso_mp(ss, constraints)
        if constraints.diso and n_diso.any():
            penalty_diso = n_diso / norm_diso_contig
        else:
            penalty_diso = np.zeros((multipanel.n_panels,))

        n_contig = calc_penalty_contig_mp(ss, constraints)
        if constraints.contig and n_contig.any():
            penalty_contig = n_contig / norm_diso_contig
        else:
            penalty_contig = np.zeros((multipanel.n_panels,))

        lampam = np.array([calc_lampam(ss[ind_panel]) \
                           for ind_panel in range(multipanel.n_panels)])

        if constraints.rule_10_percent and constraints.rule_10_Abdalla:
            penalty_10 = calc_penalty_10_ss(ss, constraints, lampam, mp=True)
        else:
            penalty_10 = calc_penalty_10_ss(ss, constraints, LPs=None)

        penalty_ipo, penalty_oopo = calc_penalty_ipo_oopo_mp(
            lampam, constraints)

        for ind_p, panel in enumerate(multipanel.panels):
            table_p.loc[ind_p, 'Penalty disorientation'] = penalty_diso[ind_p]
            table_p.loc[ind_p, 'Penalty contiguity'] = penalty_contig[ind_p]
            table_p.loc[ind_p, 'Penalty disorientation'] = penalty_diso[ind_p]
            table_p.loc[ind_p, 'Penalty contiguity'] = penalty_contig[ind_p]
            table_p.loc[ind_p, 'Penalty 10% rule'] = penalty_10[ind_p]
            table_p.loc[ind_p, 'Penalty balance'] = penalty_ipo[ind_p]
            table_p.loc[ind_p, 'Penalty out-of-plane orthotropy'] \
            = penalty_oopo[ind_p]

    if save_buckling:
        for ind_p, panel in enumerate(multipanel.panels):
            table_p.loc[ind_p, 'lambda buckling'] = buckling_factor(
                lampam=panel.lampam_target,
                mat=mat,
                n_plies=panel.n_plies,
                N_x=panel.N_x,
                N_y=panel.N_y,
                length_x=panel.length_x,
                length_y=panel.length_y,
                n_modes=10)

    append_df_to_excel(
        filename, table_p, 'Panels', index=True, header=True)
    return 0

def save_constraints_BELLA(filename, constraints):
    """
    saves the design and manufacturing constraints on Sheet [Constraints]
    """
    table_const = pd.DataFrame()
    table_const.loc[0, 'symmetry'] = constraints.sym
    table_const.loc[0, 'balance'] = constraints.bal
    table_const.loc[0, 'out-of-plane orthotropy'] = constraints.oopo
    table_const.loc[0, 'damage tolerance'] = constraints.dam_tol
    table_const.loc[0, 'dam_tol_rule'] = constraints.dam_tol_rule
    table_const.loc[0, 'covering'] = constraints.covering
    table_const.loc[0, 'n_covering'] = constraints.n_covering
    table_const.loc[0, '10% rule'] = constraints.rule_10_percent
    table_const.loc[0, '10% rule applied on LPs'] \
    = constraints.rule_10_percent and constraints.rule_10_Abdalla
    table_const.loc[0, '10% rule applied on ply percentages'] \
    = constraints.rule_10_percent and not constraints.rule_10_Abdalla
    if constraints.rule_10_percent:
        table_const.loc[0, 'percentage limit when rule applied on LPs'] \
        = constraints.percent_Abdalla * 100
        table_const.loc[0, 'percent_0'] = constraints.percent_0 * 100
        table_const.loc[0, 'percent_45'] = constraints.percent_45 * 100
        table_const.loc[0, 'percent_90'] = constraints.percent_90 * 100
        table_const.loc[0, 'percent_-45'] = constraints.percent_135 * 100
        table_const.loc[0, 'percent_+-45'] = constraints.percent_45_135 * 100
    else:
        table_const.loc[0, 'percentage limit when rule applied on LPs'] = 0
        table_const.loc[0, 'percent_0'] = 0
        table_const.loc[0, 'percent_45'] = 0
        table_const.loc[0, 'percent_90'] = 0
        table_const.loc[0, 'percent_-45'] = 0
        table_const.loc[0, 'percent_+-45'] = 0
    table_const.loc[0, 'diso'] = constraints.diso
    table_const.loc[0, 'delta_angle'] = constraints.delta_angle
    table_const.loc[0, 'contig'] = constraints.contig
    table_const.loc[0, 'n_contig'] = constraints.n_contig_c
    sets = np.array(constraints.set_of_angles, dtype=str)
    table_const.loc[0, 'fibre orientations'] = ' '.join(sets)
    table_const.loc[0, 'number fibre orientations'] \
    = constraints.n_set_of_angles
    # table_const.loc[0, 'n_plies_min'] = constraints.n_plies_min
    # table_const.loc[0, 'n_plies_max'] = constraints.n_plies_max
    table_const.loc[0, 'ply drop spacing rule'] \
    = constraints.pdl_spacing
    table_const.loc[0, 'minimum number of continuous plies between ply drops']\
    = constraints.min_drop

    table_const = table_const.transpose()

    append_df_to_excel(
        filename, table_const, 'Constraints', index=True, header=False)

def save_parameters_BELLA(filename, parameters):
    """
    saves the optimiser parameters on Sheet [Parameters]
    """
    table_param = pd.DataFrame()

    # Parameters of BELLA step 2
    table_param.loc[0, 'number of initial ply drops'] \
    = parameters.n_ini_ply_drops
    table_param.loc[0, 'minimum group size'] = parameters.group_size_min
    table_param.loc[0, 'maximum group size'] = parameters.group_size_max
    table_param.loc[0, 'time_limit_group_pdl'] = parameters.time_limit_group_pdl
    table_param.loc[0, 'time_limit_all_pdls'] = parameters.time_limit_all_pdls
    table_param.loc[0, 'global_node_limit'] \
    = parameters.global_node_limit
    table_param.loc[0, 'global_node_limit_final'] \
    = parameters.global_node_limit_final
    table_param.loc[0, 'local_node_limit'] \
    = parameters.local_node_limit
    table_param.loc[0, 'local_node_limit_final'] \
    = parameters.local_node_limit_final

    # Parameters of BELLA step 4.1
    table_param.loc[0, 'input number of plies in reference panel'] \
    = parameters.n_plies_ref_panel
    table_param.loc[0, 'repair_membrane_switch'] \
    = parameters.repair_membrane_switch
    table_param.loc[0, 'repair_flexural_switch'] \
    = parameters.repair_flexural_switch
    table_param.loc[0, 'p_A'] \
    = parameters.p_A
    table_param.loc[0, 'n_D1'] \
    = parameters.n_D1
    table_param.loc[0, 'n_D2'] \
    = parameters.n_D2
    table_param.loc[0, 'n_D3'] \
    = parameters.n_D3

    # Parameters of BELLA step 4.2
    table_param.loc[0, 'global_node_limit2'] \
    = parameters.global_node_limit2
    table_param.loc[0, 'local_node_limit2'] \
    = parameters.local_node_limit2

    # Parameters of BELLA step 4.3
    table_param.loc[0, 'global_node_limit3'] \
    = parameters.global_node_limit3
    table_param.loc[0, 'local_node_limit3'] \
    = parameters.local_node_limit3

    table_param = table_param.transpose()

    append_df_to_excel(
        filename, table_param, 'Parameters', index=True, header=False)


def save_objective_function_BELLA(filename, obj_func_param):
    """
    saves the objective function parameters on Sheet [Objective function]
    """
    table_obj_func = pd.DataFrame()

    # General parameters of BELLA
    table_obj_func.loc[0, 'optimisation problem'] = "LP matching"

#    for ind in range(12):
#        table_obj_func.loc[0, 'lampam_weightings[' + str(ind + 1) + ']'] \
#        = obj_func_param.lampam_weightings[ind]
#
#    for ind_p in range(obj_func_param.panel_weightings_ini.size):
#        table_obj_func.loc[0, 'panel_weightings_ini[' + str(ind_p + 1) + ']'] \
#        = obj_func_param.panel_weightings_ini[ind_p]

    # Penalty coefficients
    table_obj_func.loc[0, 'coeff_contig'] = obj_func_param.coeff_contig
    table_obj_func.loc[0, 'coeff_diso'] = obj_func_param.coeff_diso
    table_obj_func.loc[0, 'coeff_10'] = obj_func_param.coeff_10
    table_obj_func.loc[0, 'coeff_bal_ipo'] = obj_func_param.coeff_bal_ipo
    table_obj_func.loc[0, 'coeff_oopo'] = obj_func_param.coeff_oopo
    table_obj_func.loc[0, 'coeff_spacing'] = obj_func_param.coeff_spacing

    table_obj_func = table_obj_func.transpose()

    append_df_to_excel(filename, table_obj_func, 'Objective function',
                       index=True, header=False)