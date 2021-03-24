# -*- coding: utf-8 -*-
"""
repair for membrane properties

- repair_membrane
    repairs a laminate to improve its in-plane stiffness properties

- repair_membrane_1:
    repair for membrane properties only accounting for one panel

- repair_membrane_2:
    repair for membrane properties accounting for all the panels
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA')
from src.RELAY.repair_membrane_1_ipo import repair_membrane_1_ipo
from src.RELAY.repair_membrane_1_no_ipo import repair_membrane_1_no_ipo
from src.RELAY.repair_membrane_1_ipo_Abdalla import repair_membrane_1_ipo_Abdalla
from src.RELAY.repair_membrane_1_no_ipo_Abdalla \
import repair_membrane_1_no_ipo_Abdalla
from src.RELAY.repair_10_bal import calc_lampamA_ply_queue

def repair_membrane(
        ss, ply_queue, mini_10, in_plane_coeffs, constraints, parameters,
        obj_func_param=None, multipanel=None, lampam_target=None):
    """
    repairs a laminate to improve its in-plane stiffness properties
    """
    if not multipanel is None:
        if  not parameters.repair_membrane_switch \
        or np.isclose(np.array([0, 0, 0, 0], float), in_plane_coeffs).all():
            ss_list = [ss] # no in-plane optimisation required
            ply_queue_list = [ply_queue]
            lampamA_list = [
                calc_lampamA_ply_queue(ss, ss.size, ply_queue, constraints)]
        else:
            ss_list, ply_queue_list, lampamA_list, _ = repair_membrane_1(
                ss, ply_queue, mini_10,
                in_plane_coeffs, parameters.p_A,
                lampam_target, constraints)
        return ss_list, ply_queue_list, lampamA_list

    if not parameters.repair_membrane_switch \
    or np.isclose(np.array([0, 0, 0, 0], float), in_plane_coeffs).all():
        ss_list = [ss] # no in-plane optimisation required
        ply_queue_list = [ply_queue]
        lampamA_list = [
            calc_lampamA_ply_queue(ss, ss.size, ply_queue, constraints)]
    else:
        ss_list, ply_queue_list, lampamA_list, _ = repair_membrane_1(
            ss, ply_queue, mini_10,
            in_plane_coeffs, parameters.p_A,
            lampam_target, constraints)
    return ss_list, ply_queue_list, lampamA_list


def repair_membrane_1(
        ss, ply_queue, mini_10, in_plane_coeffs,
        p_A, lampam_target, constraints):
    """
    repair for membrane properties only accounting for one panel
    """
    if constraints.rule_10_percent and constraints.rule_10_Abdalla:
        if constraints.ipo:
            return repair_membrane_1_ipo_Abdalla(
                ss, ply_queue, in_plane_coeffs, p_A, lampam_target,
                constraints)
        return repair_membrane_1_no_ipo_Abdalla(
            ss, ply_queue, in_plane_coeffs, p_A, lampam_target, constraints)

    if constraints.ipo:
        return repair_membrane_1_ipo(
            ss, ply_queue, mini_10, in_plane_coeffs,
            p_A, lampam_target, constraints)
    return repair_membrane_1_no_ipo(
        ss, ply_queue, mini_10, in_plane_coeffs,
        p_A, lampam_target, constraints)
