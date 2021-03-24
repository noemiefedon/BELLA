# -*- coding: utf-8 -*-
"""
Optimisation of a composite laminate design
"""

__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import time
import numpy as np

sys.path.append(r'C:\BELLA')
from src.BELLA.results import BELLA_Results
#from src.BELLA.divide_panels_1 import divide_panels_1
from src.BELLA.pdl_ini import create_initial_pdls, check_pdls_ini
from src.BELLA.divide_panels import  divide_panels
from src.BELLA.ply_order import calc_ply_order
from src.BELLA.moments_of_areas import calc_mom_of_areas
from src.BELLA.lampam_matrix import calc_delta_lampams
from src.BELLA.optimiser_with_one_pdl import BELLA_optimiser_one_pdl
from src.BELLA.save_set_up import save_constraints_BELLA, save_parameters_BELLA
from src.BELLA.save_set_up import save_multipanel, save_objective_function_BELLA
from src.BELLA.save_set_up import save_materials
from src.BELLA.save_result import save_result_BELLAs
from src.divers.excel import autofit_column_widths, delete_file
from src.BELLA.format_pdl import extend_after_guide_based_blending
from src.BELLA.pdl_ini import read_pdls_excel

def BELLA_optimiser(
        multipanel, parameters, obj_func_param, constraints, filename,
        mat=None, filename_initial_pdls=None):
    """
    performs the retrieval of blended stacking sequences from
    lamination-parameter targets

    - BELLA_results: results of the optimisation

    INPUTS

    - parameters: parameters of the optimiser
    - constraints: lay-up design guidelines
    - obj_func_param: objective function parameters
    - targets: target lamination parameters and ply counts
    - mat: material properties
    - filename: name of the file where to save results
    - pdls_ini: initial ply-drop layouts (opyional)
    """
    ### initialisation
    t0 =time.time()

    delete_file(filename)

    multipanel.should_you_use_BELLA()
    multipanel.calc_weight_per_panel(mat.density_area)
    multipanel.filter_target_lampams(constraints, obj_func_param)
    multipanel.filter_lampam_weightings(constraints, obj_func_param)

    ### step 1 of BELLA: mapping the multi-panel structure to a blending strip
    print('---- Blending step 1 ----')
    multipanel.from_mp_to_blending_strip(
        constraints, parameters.n_plies_ref_panel)

    ### step 2 of BELLA: generation of initial ply drop layouts
    print('---- Blending step 2 ----')

    if filename_initial_pdls is None:
        # creation of the initial ply drop layouts
        divide_panels(multipanel, parameters, constraints)
        pdls_ini = create_initial_pdls(
            multipanel, constraints, parameters, obj_func_param)
    else:
        # read initial ply-drop layouts
        pdls_ini = read_pdls_excel(filename_initial_pdls)
        # number of initial ply drops to be tested
        parameters.n_ini_ply_drops = len(pdls_ini)
        # check the correct number of plies and panels in the ply drop layouts
        check_pdls_ini(multipanel, pdls_ini)

    ### preparation of the step 3 of BELLA

    # division of the plies of the panels into one ply group
    group_size_max = parameters.group_size_max
    parameters.group_size_max = 10000
    divide_panels(multipanel, parameters, constraints)
    parameters.group_size_max = group_size_max

    # initialisation of the results
    results = BELLA_Results(constraints, multipanel, parameters)

    # list of the orders in which plies are optimised in each panels
    ply_order = calc_ply_order(multipanel, constraints)
#    print('ply_order')
#    print(ply_order)

    # mom_areas_plus: positive ply moments of areas
    # mom_areas: ply moments of areas
    mom_areas_plus, mom_areas = calc_mom_of_areas(
        multipanel, constraints, ply_order)

    # calculation of ply partial lamination parameters
    delta_lampams = calc_delta_lampams(
        multipanel, constraints, mom_areas, ply_order)

    outer_step = 0

    while outer_step < parameters.n_ini_ply_drops:

        ### step 3 and 4 of BELLA: ply angle optimisation + laminate repair
        print('---- Blending step 3 ----')
        results_one_pdl = BELLA_optimiser_one_pdl(
            multipanel, parameters, obj_func_param, constraints, ply_order,
            mom_areas_plus, delta_lampams, pdls_ini[outer_step], mat=mat)

        results.update(outer_step, results_one_pdl)

        ## === If the stacking sequence is good enough, exit the loop
        if results_one_pdl is not None:
            if abs(results_one_pdl.obj_constraints).all() < 1e-10:
                print(f"""Low objective for the outer step {outer_step}.""")
                break

        outer_step += 1

    # To determine the best solution
    ind_mini = np.argmin(results.obj_constraints_tab)

    results.ss = results.ss_tab[ind_mini]
    results.sst = results.ss_tab_tab[ind_mini]
    results.lampam = results.lampam_tab_tab[:, ind_mini, :]
    results.n_plies_per_angle = results.n_plies_per_angle_tab[ind_mini]
    results.obj_constraints = results.obj_constraints_tab[ind_mini]
    results.obj_no_constraints = results.obj_no_constraints_tab[ind_mini]
    results.ind_mini = ind_mini
    results.time = time.time() - t0

    # save data
    save_constraints_BELLA(filename, constraints)
    if mat is not None:
        save_materials(filename, mat)
    save_parameters_BELLA(filename, parameters)
    save_objective_function_BELLA(filename, obj_func_param)
    save_multipanel(filename, multipanel, obj_func_param, mat)

    pdls = np.zeros((len(pdls_ini),
                     multipanel.n_panels,
                     multipanel.n_plies_max))
    for ind in range(len(pdls_ini)):
        pdls[ind] = np.array(
            extend_after_guide_based_blending(multipanel, pdls_ini[ind]))
    save_result_BELLAs(filename, multipanel, constraints, parameters,
                       obj_func_param, pdls, results, mat)
    autofit_column_widths(filename)

    return results
