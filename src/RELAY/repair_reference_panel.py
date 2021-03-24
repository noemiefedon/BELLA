# -*- coding: utf-8 -*-
"""
Repair strategy for reference panel
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys

sys.path.append(r'C:\BELLA')
from src.CLA.lampam_functions import calc_lampam
from src.divers.pretty_print import print_lampam, print_ss, print_list_ss
from src.BELLA.format_pdl import convert_sst_to_ss
from src.BELLA.format_pdl import convert_ss_ref_to_reduced_sst
from src.RELAY.repair_membrane import repair_membrane
from src.RELAY.repair_flexural import repair_flexural
from src.RELAY.repair_10_bal import repair_10_bal
from src.RELAY.repair_10_bal import calc_mini_10
from src.RELAY.repair_diso_contig import repair_diso_contig_list
from src.guidelines.one_stack import check_lay_up_rules


def repair_reference_panel(
        multipanel, reduced_ss, constraints, parameters, obj_func_param,
        reduced_pdl, mat=0):
    """
    repairs a reference stacking sequence to meet design and manufacturing
    guidelines

    The repair process is deterministic and attempts at conducting minimal
    modification of the original stacking sequence with a preference for
    modifying outer plies that have the least influence on out-of-plane
    properties.

    step 1: repair for the 10% rule and balance
    step 2: refinement for in-plane lamination parameter convergence
    step 3: repair for disorientation and contiguity
    step 4: refinement for out-of-plane lamination parameter convergence

    INPUTS

    - n_panels: number of panels in the entire structure
    - ss: stacking sequence of the laminate
    - multipanel: multipanel structure
    - constraints: instance of the class Constraints
    - parameters: instance of the class Parameters
    - obj_func_param: objective function parameters
    - reduced_pdl: reduced ply drop layout
    - mat: material properties
    """
#    print_list_ss(reduced_pdl, 60)

#    weight_now = mat.density_area * np.array(
#                    [outputs.ss[ind_panel].size \
#                     for ind_panel in range(multipanel.n_panels)])
#    penalty_weight_tab[outer_step] = (
#                    weight_now - weight_ref) / weight_ref

    ind_ref = multipanel.reduced.ind_ref
    ss_ref = reduced_ss[ind_ref]
    lampam_target_ref = multipanel.reduced.panels[ind_ref].lampam_target

    mini_10 = calc_mini_10(constraints, ss_ref.size)
    #--------------------------------------------------------------------------
    # step 1 / repair for the 10% rule and balance
    #--------------------------------------------------------------------------
    ss_ref, ply_queue = repair_10_bal(ss_ref, mini_10, constraints)
    #--------------------------------------------------------------------------
    # step 2 / improvement of the in-plane lamination parameter convergence
    #--------------------------------------------------------------------------
    ss_ref_list, ply_queue_list, _ = repair_membrane(
        multipanel=multipanel,
        ss=ss_ref,
        ply_queue=ply_queue,
        mini_10=mini_10,
        in_plane_coeffs=multipanel.reduced.panels[ind_ref].lampam_weightingsA,
        parameters=parameters,
        obj_func_param=obj_func_param,
        constraints=constraints,
        lampam_target=lampam_target_ref)
    #--------------------------------------------------------------------------
    # step 3 / repair for disorientation and contiguity
    #--------------------------------------------------------------------------
    ss_ref, completed_inward, completed_outward, ind = repair_diso_contig_list(
        ss_ref_list, ply_queue_list, constraints,
        parameters.n_D1)
    if not completed_outward:

        reduced_sst = convert_ss_ref_to_reduced_sst(
                ss_ref, reduced_pdl=reduced_pdl,
                ind_ref=multipanel.reduced.ind_ref,
                reduced_ss_before=reduced_ss)

        reduced_lampam = calc_lampam(reduced_ss, constraints)
        return False, reduced_lampam, reduced_sst, reduced_ss
    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------
    reduced_sst = convert_ss_ref_to_reduced_sst(
        ss_ref, reduced_pdl=reduced_pdl, ind_ref=ind_ref,
        reduced_ss_before=reduced_ss)
    reduced_ss = convert_sst_to_ss(reduced_sst)
    reduced_lampam = calc_lampam(reduced_ss, constraints)
    #--------------------------------------------------------------------------
    # step 4 / improvement of the out-of-plane lamination parameter convergence
    #--------------------------------------------------------------------------
    ss_ref = repair_flexural(
        ss=ss_ref,
        lampam_target=lampam_target_ref,
        out_of_plane_coeffs=multipanel.reduced.panels[
            ind_ref].lampam_weightingsD,
        parameters=parameters,
        constraints=constraints,
        multipanel=multipanel)

    reduced_sst = convert_ss_ref_to_reduced_sst(
                ss_ref, reduced_pdl=reduced_pdl,
                ind_ref=multipanel.reduced.ind_ref,
                reduced_ss_before=reduced_ss)

    reduced_ss = convert_sst_to_ss(reduced_sst)
    reduced_lampam = calc_lampam(reduced_ss, constraints)

    check_lay_up_rules(ss_ref, constraints)

    if (reduced_ss[ind_ref] != ss_ref).any():
        print(ss_ref - reduced_ss[ind_ref])
        raise Exception("""
Reference stacking sequence in reduced_ss different from ss_ref""")

    return True, reduced_lampam, reduced_sst, reduced_ss
