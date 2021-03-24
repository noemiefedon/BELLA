# -*- coding: utf-8 -*-
"""
Repair strategy for guide based blending
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA')
#from src.divers.pretty_print import print_lampam, print_ss, print_list_ss
from src.BELLA.format_pdl import convert_ss_to_sst
from src.RELAY.thick_to_thin import repair_thick_to_thin
from src.RELAY.thin_to_thick import repair_thin_to_thick
from src.RELAY.repair_reference_panel import repair_reference_panel
from src.guidelines.one_stack import check_ply_drop_rules

def repair_mp(multipanel, reduced_ss, constraints, parameters, obj_func_param,
              reduced_pdl, mat=0):
    """
    repairs a multi-panel design to meet design and manufacturing guidelines
    and evaluates the performance of the repaired stacking sequence.

    The repair process is deterministic and attempts at conducting minimal
    modification of the original stacking sequence with a preference for
    modifying outer plies that have the least influence on out-of-plane
    properties.

    step 1:
        repair of the reference panel stacking sequence

    step 2:
        repair of the other panels by re-designing the ply drop layout

    INPUTS

    - n_panels: number of panels in the entire structure
    - reduced_ss: stacking sequence of the laminate
    - multipanel: multipanel structure
    - constraints: instance of the class Constraints
    - parameters: instance of the class Parameters
    - obj_func_param: objective function parameters
    - reduced_pdl: ply drop layout
    - mat: material properties
    """
    #--------------------------------------------------------------------------
    # step 1 / reference panel repair
    #--------------------------------------------------------------------------
    print('---- Blending step 4.1 ----')
    success, reduced_lampam, reduced_sst, reduced_ss = repair_reference_panel(
        multipanel, reduced_ss, constraints, parameters, obj_func_param,
        reduced_pdl, mat=0)
    if not success:
        print('Blending step 4.1 unsuccessful')
        reduced_sst = convert_ss_to_sst(reduced_ss, reduced_pdl)
        return False, reduced_ss, reduced_sst, reduced_pdl, 1

    #--------------------------------------------------------------------------
    # step 2 / re-optimise the ply drop layout - thick-to-thin repair
    #--------------------------------------------------------------------------
    print('---- Blending step 4.2 ----')
#    print('reduced_sst', reduced_sst.shape)
#    print_list_ss(reduced_sst[:,:reduced_sst.shape[1] // 2])
#    print('SS_ref')
#    print_ss(ss_ref[:ss_ref.size // 2])
    ss_ref = np.copy(reduced_ss[multipanel.reduced.ind_ref])
    success, reduced_sst, reduced_lampam, reduced_ss = repair_thick_to_thin(
        reduced_lampam, reduced_sst, reduced_ss, multipanel,
        parameters, obj_func_param, constraints, mat=mat)
    if not success:
        print('Blending step 4.2 unsuccessful')
        reduced_sst = convert_ss_to_sst(reduced_ss, reduced_pdl)
        return False, reduced_ss, reduced_sst, reduced_pdl, 2
    # check that the reference panel stacking sequence has not been changed
    if (reduced_ss[multipanel.reduced.ind_ref] != ss_ref).any():
        raise Exception("""
Reference stacking sequence modified during blending step 4.2""")

    #--------------------------------------------------------------------------
    # step 3 / re-optimise the ply drop layout - thin-to-thck repair
    #--------------------------------------------------------------------------
    print('---- Blending step 4.3 ----')
    success, reduced_sst, reduced_lampam, reduced_ss = repair_thin_to_thick(
        reduced_lampam, reduced_sst, reduced_ss, multipanel,
        parameters, obj_func_param, constraints, mat=mat)
    if not success:
        print('Blending step 4.3 unsuccessful')
        reduced_sst = convert_ss_to_sst(reduced_ss, reduced_pdl)
        return False, reduced_ss, reduced_sst, reduced_pdl, 3

    # check that the reference panel stacking sequence has not been changed
    if (reduced_ss[multipanel.reduced.ind_ref] != ss_ref).any():
        raise Exception("""
Reference stacking sequence modified during blending step 4.3""")

    #--------------------------------------------------------------------------
    # return result
    #--------------------------------------------------------------------------
#    check_ply_drop_rules(reduced_sst, multipanel, constraints)

    # test for the ply counts
    for ind_panel, panel in enumerate(multipanel.reduced.panels):
        if reduced_ss[ind_panel].size != panel.n_plies:
            raise Exception("""Wrong ply counts in the laminate.""")

#    print('Blending step 4 successful')
    return True, reduced_ss, reduced_sst, reduced_pdl, 4
