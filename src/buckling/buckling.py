# -*- coding: utf-8 -*-
"""
Functions to calculate buckling resistance

- buckling_margin
    calculates buckling margins from src.buckling factors

- buckling_factor_ss
    calculates the critical buckling factor of a simply-supported orthotropic
    laminate plate based on the stacking sequqnce

- buckling_factor_lampam
    calculates the critical buckling factor of a simply-supported orthotropic
    laminate plate based on lamination parameters

- buckling_factor_m_n
    returns the buckling factor of a simply-supported orthotropic laminate
    for a specific buckling mode, based on out-of-plane stiffnesses
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import math
import numpy as np
from src.BELLA.materials import Material
from src.CLA.lampam_functions import calc_lampam

def buckling_margin(buckling_factor):
    'calculates buckling margins from src.buckling factors'
    return 100*(buckling_factor - 1)


def buckling_factor_ss(ss, N_x, N_y, length_x, length_y, mat, n_modes=10):
    """
    calculates the critical buckling factor of a simply-supported orthotropic
    laminate plate based on the stacking sequqnce

    INPUT

    - ss: stacking sequence of the laminate
    - mat: material properties
    - length_x: plate dimensions (x-direction)
    - length_y: plate dimensions (y-direction)
    - N_x: compressive loading intensity in the x-direction
    - N_y: compressive loading intensity in the x-direction
    - n_modes: number of buckling modes to be tested
    """
    lampam_1 = calc_lampam(ss)
    n_plies = ss.size
    return buckling_factor(
        lampam_1, mat, n_plies, N_x=N_x, N_y=N_y,
        length_x=length_x, length_y=length_y, n_modes=n_modes)


def buckling_factor(
        lampam, mat, n_plies, N_x=0, N_y=0, length_x=1, length_y=1, n_modes=10):
    """
    calculates the critical buckling factor of a simply-supported orthotropic
    laminate plate based on lamination parameters

    INPUT

    - lampam: lamination parameters
    - mat: material properties
    - length_x: plate dimensions (x-direction)
    - length_y: plate dimensions (y-direction)
    - N_x: compressive loading intensity in the x-direction
    - N_y: compressive loading intensity in the x-direction
    - n_modes: number of buckling modes to be tested
    """
    buck = np.zeros((n_modes, n_modes), dtype=float)
    a = (1/12)*((n_plies*mat.ply_t)**3)
    if lampam.size == 12:
        D11 = a * (mat.U1 + mat.U2*lampam[8] + mat.U3*lampam[9])
        D12 = a *(-mat.U3*lampam[9] + mat.U4)
        D22 = a *(mat.U1 - mat.U2*lampam[8] + mat.U3*lampam[9])
        D66 = a *(-mat.U3*lampam[9] + mat.U5)
    else:
        D11 = a * (mat.U1 + mat.U2*lampam[0] + mat.U3*lampam[1])
        D12 = a *(-mat.U3*lampam[1] + mat.U4)
        D22 = a *(mat.U1 - mat.U2*lampam[0] + mat.U3*lampam[1])
        D66 = a *(-mat.U3*lampam[0] + mat.U5)
    for mode_m in range(n_modes):
        for mode_n in range(n_modes):
            buck[mode_m, mode_n] = buckling_factor_m_n(
                D11, D12, D22, D66,
                mode_m=mode_m + 1,
                mode_n=mode_n + 1,
                N_x=N_x, N_y=N_y,
                length_x=length_x,
                length_y=length_y)
    return np.min(buck)


def buckling_factor_m_n(
        D11, D12, D22, D66, mode_m=1, mode_n=1,
        N_x=0, N_y=0, length_x=1, length_y=1):
    """
    returns the buckling factor of a simply-supported orthotropic laminate
    for a specific buckling mode, based on out-of-plane stiffnesses

    INPUT

    - D11, D12, D22 and D66: out-of-plane stifnesses
    - length_x: plate dimensions (x-direction)
    - length_y: plate dimensions (y-direction)
    - N_x: compressive loading intensity in the x-direction
    - N_y: compressive loading intensity in the x-direction
    - mode_m: buckling mode (number of half-waves) in the x direction
    - mode_n: buckling mode (number of half-waves) in the y direction
    """
    alpha = (mode_m/length_x)**2
    beta = (mode_n/length_y)**2
    return (math.pi**2)*(D11*alpha**2 \
        + 2*(D12 + 2*D66)*alpha*beta \
        + D22*beta**2)/(alpha*abs(N_x) + beta*abs(N_y))


if __name__ == "__main__":
    print('*** Test for the functions buckling_factor_ss ***\n')
    # Elastic modulus in the fibre direction in Pa
    E11 = 20.5/1.45038e-10  # 141 GPa
    # Elastic modulus in the transverse direction in Pa
    E22 = 1.31/1.45038e-10  # 9.03 GPa
    # Poisson's ratio relating transverse deformation and axial loading (-)
    nu12 = 0.32
    # In-plane shear modulus in Pa
    G12 = 0.62/1.45038e-10 # 4.27 GPa
    # Density in g/m2
    density_area = 300.5
    # Ply thickness in m
    ply_t = (25.40/1000)*0.0075 # 0.191 mmm
    mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12,
                   density_area=density_area, ply_t=ply_t)

    # panel 17
    ss = np.array([45, 60, 45, 90, 90, -45, -60, -45, 0]) # panel 17 Irisarri 1
    ss = np.array([60, -60, 60, -60, 60, -60, -75, 60, 75]) # panel 17 Adams
    ss = np.array([45, -45, 45, -45, 45, -45, 60, -60, 60, -30]) # panel 17 Serestra
    N_x = 1000*0.175127*320
    N_y = 1000*0.175127*180
    length_x = 20*(25.40/1000)
    length_y = 12*(25.40/1000)
#    # panel 7
#    N_x = 1000*0.175127*290
#    N_y = 1000*0.175127*195
#    length_x = 20*(25.40/1000)
#    length_y = 12*(25.40/1000)
#    ss = np.array([45, 60, 45, 90, 90, -45, -60, -45, 0]) # panel 7 Irisarri 1
#    ss = np.array([60, -60, 60, -60, 60, -60, -75, 60, 75]) # panel 7 Adams
#    ss = np.array([45, -45, 45, -45, 45, -45, 60, -60, 60, -30]) # panel 7 Serestra
    # panel 5
    ss = np.array([45, -45, 45, -45, 45, -45, 60, -60]) # panel 5 Serestra
    ss = np.array([-60, -60, -60, -60, -60, -60, -60, -60]) # my test
    N_x = 1000*0.175127*210
    N_y = 1000*0.175127*100
    length_x = 20*(25.40/1000)
    length_y = 12*(25.40/1000)
    ss = np.hstack((ss, np.flip(ss, axis=0)))
    print(ss)
    buck = buckling_factor_ss(ss, N_x, N_y, length_x, length_y, mat, n_modes=10)
    print(f'Buckling factor : {buck }\n')
    print(f'Buckling margin : {buckling_margin(buck)}\n')


    print('*** Test for the functions buckling_factor ***\n')
    # Lamination paraneters calculation
    ss = np.array([-30, 30, 45, -45, 45, 60, -60,
                   60, -60, 60, -60, -75, 60, 75])
    n_plies = ss.size
    lampam = calc_lampam(ss)
    # Elastic modulus in the fibre direction in Pa
    E11 = 20.5/1.45038e-10  # 141 GPa
    # Elastic modulus in the transverse direction in Pa
    E22 = 1.31/1.45038e-10  # 9.03 GPa
    # Poisson's ratio relating transverse deformation and axial loading (-)
    nu12 = 0.32
    # In-plane shear modulus in Pa
    G12 = 0.62/1.45038e-10 # 4.27 GPa
    # Density in g/m2
    density_area = 300.5
    # Ply thickness in m
    ply_t = (25.40/1000)*0.0075 # 0.191 mmm
    mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12,
                   density_area=density_area, ply_t=ply_t)
    # Loading intensities
    N_x = 1000*0.175127*375
    N_y = 1000*0.175127*360
    length_x = 18*(25.40/1000)
    length_y = 24*(25.40/1000)
    buck = buckling_factor(lampam, mat, n_plies, N_x=N_x, N_y=N_y,
                           length_x=length_x, length_y=length_y)
    print(f'Buckling factor : {buck}\n')
    print(f'Buckling margin : {buckling_margin(buck)}\n')

    print('*** Test for the functions buckling_factor_m_n ***\n')
    # Lamination paraneters calculation
    ss = np.array([-30, 30, 45, -45, 45, 60, -60,
                   60, -60, 60, -60, -75, 60, 75])
    n_plies = ss.size
    lampam = calc_lampam(ss)
    # Elastic modulus in the fibre direction in Pa
    E11 = 20.5/1.45038e-10  # 141 GPa
    # Elastic modulus in the transverse direction in Pa
    E22 = 1.31/1.45038e-10  # 9.03 GPa
    # Poisson's ratio relating transverse deformation and axial loading (-)
    nu12 = 0.32
    # In-plane shear modulus in Pa
    G12 = 0.62/1.45038e-10 # 4.27 GPa
    # Density in g/m2
    density_area = 300.5
    # Ply thickness in m
    ply_t = (25.40/1000)*0.0075 # 0.191 mmm
    mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12,
                   density_area=density_area, ply_t=ply_t)
    # Loading intensities
    N_x = 1000*0.175127*375
    N_y = 1000*0.175127*360
    length_x = 18*(25.40/1000)
    length_y = 24*(25.40/1000)
    # Buckling modes
    mode_m = 1
    mode_n = 1
    # out-of-plane stiffnesses
    a = (1/12)*((n_plies*mat.ply_t)**3)
    D11 = a * (mat.U1 + mat.U2*lampam[8] + mat.U3*lampam[9])
    D12 = a *(-mat.U3*lampam[9] + mat.U4)
    D22 = a *(mat.U1 - mat.U2*lampam[8] + mat.U3*lampam[9])
    D66 = a *(-mat.U3*lampam[9] + mat.U5)
    buck = buckling_factor_m_n(D11, D12, D22, D66,
                               mode_m=mode_m, mode_n=mode_n, N_x=N_x, N_y=N_y,
                               length_x=length_x, length_y=length_y)
    print(f'Buckling factor : {buck}\n')
    print(f'Buckling margin : {buckling_margin(buck)}\n')


    print('*** Test for the functions buckling_factor ***\n')
    # panel 18 terence
    n_plies = 22
    lampam = np.array([0, 0, 0, 0,    0, 0, 0, 0,    -0.469, -0.335, 0, 0])
    # Elastic modulus in the fibre direction in Pa
    E11 = 20.5/1.45038e-10  # 141 GPa
    # Elastic modulus in the transverse direction in Pa
    E22 = 1.31/1.45038e-10  # 9.03 GPa
    # Poisson's ratio relating transverse deformation and axial loading (-)
    nu12 = 0.32
    # In-plane shear modulus in Pa
    G12 = 0.62/1.45038e-10 # 4.27 GPa
    # Density in g/m2
    density_area = 300.5
    # Ply thickness in m
    ply_t = (25.40/1000)*0.0075 # 0.191 mmm
    mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12,
                   density_area=density_area, ply_t=ply_t)
    # Loading intensities
    N_x = 1000*0.175127*(-300)
    N_y = 1000*0.175127*(-410)
    length_x = 20*(25.40/1000)
    length_y = 12*(25.40/1000)
    buck = buckling_factor(lampam, mat, n_plies, N_x=N_x, N_y=N_y,
                           length_x=length_x, length_y=length_y)
    print(f'1 - buckling factor : {1 - buck}\n')

    # panel 18 moi
    n_plies = 22
    lampam = np.array([0, 0, 0, 0,    0, 0, 0, 0,    -0.416, -0.451, 0, 0])
    # Elastic modulus in the fibre direction in Pa
    E11 = 20.5/1.45038e-10  # 141 GPa
    # Elastic modulus in the transverse direction in Pa
    E22 = 1.31/1.45038e-10  # 9.03 GPa
    # Poisson's ratio relating transverse deformation and axial loading (-)
    nu12 = 0.32
    # In-plane shear modulus in Pa
    G12 = 0.62/1.45038e-10 # 4.27 GPa
    # Density in g/m2
    density_area = 300.5
    # Ply thickness in m
    ply_t = (25.40/1000)*0.0075 # 0.191 mmm
    mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12,
                   density_area=density_area, ply_t=ply_t)
    # Loading intensities
    N_x = 1000*0.175127*(-300)
    N_y = 1000*0.175127*(-410)
    length_x = 20*(25.40/1000)
    length_y = 12*(25.40/1000)
    buck = buckling_factor(lampam, mat, n_plies, N_x=N_x, N_y=N_y,
                           length_x=length_x, length_y=length_y)
    print(f'1 - buckling factor : {1 - buck}\n')