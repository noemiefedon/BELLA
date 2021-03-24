# -*- coding: utf-8 -*-
"""
Functions related to orthotropy requirements

- calc_penalty_ipo
    calculates penalties for the balance constraint based on lamination
    parameters

- ipo_param_1_12
    calculates the twelve laminate in-plane orthotropy parameters

- ipo_param_1_6
    calculates the first six laminate in-plane orthotropy parameters

- ipo_param_7_12
    calculates the last six laminate in-plane orthotropy parameters

- calc_penalty_ipo_param
    calculates penalties for in-plane orthotropy based on in-plane orthotropy
    parameters

- calc_penalty_ipo_oopo_mp
    calculates penalties for in-plane orthotropy based on in-plane orthotropy
    lamination parameters for a multi-panel structure

- calc_penalty_ipo_oopo_ss
    calculates penalties for in-plane and out-of plane orthotropy based
    lamination parameters for a single-panel structure

- calc_penalty_oopo_ss
    calculates penalties for out-of plane orthotropy based lamination
    parameters for a single-panel structure
    """
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np
sys.path.append(r'C:\BELLA')
from src.CLA.lampam_functions import calc_lampam
from src.BELLA.materials import Material
from src.BELLA.parameters import Parameters
from src.BELLA.constraints import Constraints
from src.BELLA.panels import Panel
from src.BELLA.multipanels import MultiPanel
from src.BELLA.divide_panels import divide_panels

def ipo_param_1_12(lampam, material, sym):
    'calculates the twelve laminate in-plane orthotropy parameters'
    # extensional stifness matrix
    A11 = (material.U1 + material.U2*lampam[0] + material.U3*lampam[1])
    A12 = (- material.U3*lampam[1] + material.U4)
    A22 = (material.U1 - material.U2*lampam[0] + material.U3*lampam[1])
    A66 = (- material.U3*lampam[1] + material.U5)
    A16 = (0.5*material.U2*lampam[2] + material.U3*lampam[3])
    A26 = (0.5*material.U2*lampam[2] - material.U3*lampam[3])
    param1 = A16/A11
    param2 = A16/A12
    param3 = A16/A66
    param4 = A26/A12
    param5 = A26/A22
    param6 = A26/A66
    if sym:
        param7 = (A12*A26 - A22*A16)/(A22*A66 - A26*A26)
        param8 = (A12*A26 - A22*A16)/(A16*A26 - A12*A66)
        param9 = (A12*A26 - A22*A16)/(A11*A22 - A12*A12)
        param10 = (A12*A16 - A11*A26)/(A16*A26 - A12*A66)
        param11 = (A12*A16 - A11*A26)/(A11*A66 - A16*A16)
        param12 = (A12*A16 - A11*A26)/(A11*A22 - A12*A12)
    else:
        # coupling stifness matrix
        B11 = (0.25)*(material.U2*lampam[4] + material.U3*lampam[5])
        B12 = (0.25)*(- material.U3*lampam[5])
        B22 = (0.25)*(- material.U2*lampam[4] + material.U3*lampam[5])
        B66 = (0.25)*(- material.U3*lampam[5])
        B16 = (0.25)*(0.5*material.U2*lampam[6] + material.U3*lampam[7])
        B26 = (0.25)*(0.5*material.U2*lampam[6] - material.U3*lampam[7])
        # bend/twist stifness matrix
        D11 = (1/12)*(material.U1 + material.U2*lampam[8] + material.U3*lampam[9])
        D12 = (1/12)*(- material.U3*lampam[9] + material.U4)
        D22 = (1/12)*(material.U1 - material.U2*lampam[8] + material.U3*lampam[9])
        D66 = (1/12)*(- material.U3*lampam[9] + material.U5)
        D16 = (1/12)*(0.5*material.U2*lampam[10] + material.U3*lampam[11])
        D26 = (1/12)*(0.5*material.U2*lampam[10] - material.U3*lampam[11])
        A = np.array([[A11, A12, A16],
                      [A12, A22, A26],
                      [A16, A26, A66]])
        B = np.array([[B11, B12, B16],
                      [B12, B22, B26],
                      [B16, B26, B66]])
        D = np.array([[D11, D12, D16],
                      [D12, D22, D26],
                      [D16, D26, D66]])
        # reduced menbrane compliance
        a = np.linalg.inv(A - B@(np.linalg.inv(D))@B)
        a11 = a[0, 0]
        a12 = a[0, 1]
        a22 = a[1, 1]
        a66 = a[2, 2]
        a16 = a[0, 2]
        a26 = a[1, 2]
        param7 = a16/a11
        param8 = a16/a12
        param9 = a16/a66
        param10 = a26/a12
        param11 = a26/a22
        param12 = a26/a66
    return abs(np.array([param1, param2, param3,
                         param4, param5, param6,
                         param7, param8, param9,
                         param10, param11, param12]))

def ipo_param_1_6(lampam, material, sym):
    'calculates the first six laminate in-plane orthotropy parameters'
    # extensional stifness matrix
    A11 = (material.U1 + material.U2*lampam[0] + material.U3*lampam[1])
    A12 = (- material.U3*lampam[1] + material.U4)
    A22 = (material.U1 - material.U2*lampam[0] + material.U3*lampam[1])
    A66 = (- material.U3*lampam[1]+ material.U5)
    A16 = (0.5*material.U2*lampam[2] + material.U3*lampam[3])
    A26 = (0.5*material.U2*lampam[2] - material.U3*lampam[3])
    param1 = A16/A11
    param2 = A16/A12
    param3 = A16/A66
    param4 = A26/A12
    param5 = A26/A22
    param6 = A26/A66
    param1 = A16/A11
    param2 = A16/A12
    param3 = A16/A66
    param4 = A26/A12
    param5 = A26/A22
    param6 = A26/A66
    return abs(np.array([param1, param2, param3,
                         param4, param5, param6]))


def ipo_param_7_12(lampam, material, sym):
    'calculates the last six laminate in-plane orthotropy parameters'
    # extensional stifness matrix
    A11 = (material.U1 + material.U2*lampam[0] + material.U3*lampam[1])
    A12 = (- material.U3*lampam[1] + material.U4)
    A22 = (material.U1 - material.U2*lampam[0] + material.U3*lampam[1])
    A66 = (- material.U3*lampam[1] + material.U5)
    A16 = (0.5*material.U2*lampam[2] + material.U3*lampam[3])
    A26 = (0.5*material.U2*lampam[2] - material.U3*lampam[3])
    if sym:
        param7 = (A12*A26 - A22*A16)/(A22*A66 - A26*A26)
        param8 = (A12*A26 - A22*A16)/(A16*A26 - A12*A66)
        param9 = (A12*A26 - A22*A16)/(A11*A22 - A12*A12)
        param10 = (A12*A16 - A11*A26)/(A16*A26 - A12*A66)
        param11 = (A12*A16 - A11*A26)/(A11*A66 - A16*A16)
        param12 = (A12*A16 - A11*A26)/(A11*A22 - A12*A12)
    else:
        # coupling stifness matrix
        B11 = (0.25)*(material.U2*lampam[4] + material.U3*lampam[5])
        B12 = (0.25)*(- material.U3*lampam[5])
        B22 = (0.25)*(- material.U2*lampam[4] + material.U3*lampam[5])
        B66 = (0.25)*(- material.U3*lampam[5])
        B16 = (0.25)*(0.5*material.U2*lampam[6] + material.U3*lampam[7])
        B26 = (0.25)*(0.5*material.U2*lampam[6] - material.U3*lampam[7])
        # bend/twist stifness matrix
        D11 = (1/12)*(material.U1 + material.U2*lampam[8] + material.U3*lampam[9])
        D12 = (1/12)*(- material.U3*lampam[9] + material.U4)
        D22 = (1/12)*(material.U1 - material.U2*lampam[8] + material.U3*lampam[9])
        D66 = (1/12)*(- material.U3*lampam[9] + material.U5)
        D16 = (1/12)*(0.5*material.U2*lampam[10] + material.U3*lampam[11])
        D26 = (1/12)*(0.5*material.U2*lampam[10] - material.U3*lampam[11])
        A = np.array([[A11, A12, A16],
                      [A12, A22, A26],
                      [A16, A26, A66]])
        B = np.array([[B11, B12, B16],
                      [B12, B22, B26],
                      [B16, B26, B66]])
        D = np.array([[D11, D12, D16],
                      [D12, D22, D26],
                      [D16, D26, D66]])
        # reduced menbrane compliance
        a = np.linalg.inv(A - B@(np.linalg.inv(D))@B)
        a11 = a[0, 0]
        a12 = a[0, 1]
        a22 = a[1, 1]
        a66 = a[2, 2]
        a16 = a[0, 2]
        a26 = a[1, 2]
        param7 = a16/a11
        param8 = a16/a12
        param9 = a16/a66
        param10 = a26/a12
        param11 = a26/a22
        param12 = a26/a66
    return abs(np.array([param7, param8, param9,
                         param10, param11, param12]))


def calc_penalty_ipo_param(param, threshold):
    """
    calculates penalties for in-plane orthotropy based on in-plane orthotropy
    parameters
    """
    return max(np.array([max(0, abs(param[ii]))
                         for ii in range(param.size)]))/param.size
#    return sum(np.array([max(0, abs(param[ii]) - threshold)/threshold\
#                         for ii in range(param.size)]))/param.size


def calc_penalty_ipo(lampam, cummul_areas=1):
    """
    calculates penalties for the balance constraint based on lamination
    parameters
    """
    if lampam.ndim == 2:
        n_panels = lampam.shape[0]
        penalties_ipo = np.zeros((n_panels,), dtype=float)
        for ind_panel in range(n_panels):
            penalties_ipo[ind_panel] = (
                abs(lampam[ind_panel][2]) + abs(lampam[ind_panel][3])) / 2
    else:
        penalties_ipo = (abs(lampam[2]) + abs(lampam[3])) / 2
    return cummul_areas * penalties_ipo


def calc_penalty_ipo_oopo_mp(
        lampam,
        constraints,
        penalty_ipo_switch=True,
        parameters=None,
        mat=0,
        cummul_areas=1,
        cummul_sec_mom_areas=1):
    """
    calculates penalties for in-plane orthotropy based on in-plane orthotropy
    lamination parameters for a multi-panel structure

    INPUTS

    - lampam: panel lamination parameters
    - mat: material properties of the laminae
    - constraints: design and maufacturing constraints
    - parameters: optimiser parameters
    - cummul_areas: sum of the areas of the plies retrieved so far + the
    current plies
    - cummul_sec_mom_areas: sum of the second moments of areas of the plies
    retrieved so far + the current plies
    """
    if parameters is None:
        calculate_penalty = True

    n_panels = lampam.shape[0]
    penalties_ipo = np.zeros((n_panels,), dtype=float)
    penalties_oopo = np.zeros((n_panels,), dtype=float)
    if calculate_penalty:
        for ind_panel in range(n_panels):
            # penalty for in-plane orthotropy
            if constraints.ipo and penalty_ipo_switch:
                penalties_ipo[ind_panel] = (
                    abs(lampam[ind_panel][2]) + abs(lampam[ind_panel][3])) / 2
            # penalty for out-of-plane orthotropy
            if constraints.oopo:
                penalties_oopo[ind_panel] = (
                    abs(lampam[ind_panel][10]) + abs(lampam[ind_panel][11])) / 2
    return cummul_areas * penalties_ipo, cummul_sec_mom_areas * penalties_oopo



def calc_penalty_oopo_ss(lampam, constraints, cummul_sec_mom_areas=1):
    """
    calculates penalties for out-of plane orthotropy based lamination
    parameters for a single-panel structure

    INPUTS

    - lampam: lamination parameters
    - constraints: design and maufacturing constraints
    - cummul_sec_mom_areas: sum of the second moments of areas of the plies
    retrieved so far + the current plies
    """
        
    if (isinstance(lampam, list) and len(lampam) > 1) or lampam.ndim == 2:
        if (isinstance(lampam, list) and len(lampam) > 1):
            n_ss = len(lampam)
        else:
            n_ss = lampam.shape[0]
            
        penalties_oopo = np.zeros((n_ss,), dtype=float)
        if constraints.oopo:
            for ind_ss in range(n_ss):
                penalties_oopo[ind_ss] = (
                    abs(lampam[ind_ss][10]) + abs(lampam[ind_ss][11])) / 2
    else:
        penalties_oopo = 0
        if constraints.oopo:
            penalties_oopo = (abs(lampam[10]) + abs(lampam[11])) / 2

    return cummul_sec_mom_areas * penalties_oopo



if __name__ == "__main__":
    ss = np.array([0, 45, 45, -45, -45, 90, 45, 90])
    ss = np.hstack((ss, np.flip(ss, axis=0)))
    lampam = calc_lampam(ss)
    E11 = 130e9
    E22 = 9e9
    nu12 = 0.3
    G12 = 4e9
    threshold = 0.01
    mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12)
    sym = False
    print("""*** Test for the functions calc_penalty_ipo_param and ipo_param***\n""")
    param = ipo_param_1_6(lampam, mat, sym)
    print(f'In-plane orthotropy parameters = \n {param[0:6]}\n')
    print(f'calc_penalty_ipo : {calc_penalty_ipo_param(param, threshold)}\n')
    param = ipo_param_7_12(lampam, mat, sym)
    print(f'In-plane orthotropy parameters = \n {param[0:6]}\n')
    print(f'calc_penalty_ipo : {calc_penalty_ipo_param(param, threshold)}\n')
    param = ipo_param_1_12(lampam, mat, sym)
    print(f'In-plane orthotropy parameters = \n {param[0:6]} \n{param[6:12]}\n')
    print(f'calc_penalty_ipo : {calc_penalty_ipo_param(param, threshold)}\n')


    print('\n*** Test for the functions calc_penalty_ipo_oopo_mp***\n')
    constraints = Constraints(bal=True, oopo=True)
    parameters = Parameters(constraints)
    lampam = 0.22*np.ones((2, 12))
    lampam_target = 0.11*np.ones((2, 12))
    lampam_weightings = np.ones((12,), dtype=float)
    E11 = 130e9
    E22 = 9e9
    nu12 = 0.3
    G12 = 4e9
    mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12)
    panel_1 = Panel(ID=1,
                    neighbour_panels=[2],
                    lampam_target=0.66*np.ones((12,), dtype=float),
                    n_plies=12,
                    area=1,
                    constraints=constraints)
    panel_2 = Panel(ID=2,
                    lampam_target=0.66*np.ones((12,), dtype=float),
                    n_plies=10,
                    area=1,
                    constraints=constraints)
    panel_1.lampam_weightings = np.ones((12,), dtype=float)
    panel_2.lampam_weightings = np.ones((12,), dtype=float)
    boundaries = np.array([[0, 1]])
    multipanel = MultiPanel(panels=[panel_1, panel_2])
    print(calc_penalty_ipo_oopo_mp(
        lampam,
        constraints=constraints,
        parameters=parameters,
        mat=mat,
        cummul_areas=4,
        cummul_sec_mom_areas=20))


