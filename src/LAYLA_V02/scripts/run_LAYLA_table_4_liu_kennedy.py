# -*- coding: utf-8 -*-
"""
LAYLA retrieves the LAminate LAY-ups from lamination parameters corresponding
to the Table 4 in the publication:

    "Two-level layup optimization of composite laminate using lamination
    parameters " Liu X Featherston C Kennedy D
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np
import numpy.matlib
import pandas as pd
import time
import sys
sys.path.append(r'C:\BELLA_and_LAYLA')
from src.LAYLA_V02.targets import Targets
from src.LAYLA_V02.parameters import Parameters
from src.LAYLA_V02.constraints import Constraints
from src.LAYLA_V02.optimiser import LAYLA_optimiser
from src.LAYLA_V02.objectives import objectives

from src.BELLA.materials import Material

from src.CLA.lampam_functions import calc_lampam
from src.CLA.ABD import A_from_lampam, B_from_lampam, D_from_lampam

from src.guidelines.one_stack import check_lay_up_rules

from src.divers.excel import autofit_column_widths
from src.divers.excel import delete_file, append_df_to_excel
from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

from src.LAYLA_V02.save_set_up import save_constraints_LAYLA
from src.LAYLA_V02.save_set_up import save_parameters_LAYLA_V02
from src.LAYLA_V02.save_set_up import save_materials


result_filename = 'LAYLA_vs_global_layerwise_optimisation.xlsx'
delete_file(result_filename)

#==============================================================================
# design and manufacturing constraints
#==============================================================================
set_of_angles = np.array([-45, 0, 45, 90], dtype=int)
#set_of_angles = np.array([-45, 0, 45, 90, +30, -30, +60, -60], dtype=int)

# rule 1: one outer ply at + or -45 deg at laminate surfaces
# rule 2: [+45, -45] or [-45, +45] plies at laminate surfaces
# rule 3: [+45, -45], [+45, +45], [-45, -45] or [-45, +45] plies at laminate
dam_tol_rule = 3

combine_45_135 = False
percent_0 = 10 # percentage used in the 10% rule for 0 deg plies
percent_45 =  10 # percentage used in the 10% rule for +45 deg plies
percent_90 = 10 # percentage used in the 10% rule for 90 deg plies
percent_135 = 10 # percentage used in the 10% rule for -45 deg plies
percent_45_135 = 10 # percentage used in the 10% rule for +-45 deg plies

delta_angle = 45

n_contig = 4

#==============================================================================
# Material properties
#==============================================================================
# Elastic modulus in the fibre direction (Pa)
E11 = 128e9
# Elastic modulus in the transverse direction (Pa)
E22 = 10.3e9
# Poisson's ratio relating transverse deformation and axial loading (-)
nu12 = 0.3
# In-plane shear modulus (Pa)
G12 = 6e9
mat_prop = Material(E11 = E11, E22 = E22, G12 = G12, nu12 = nu12)

#==============================================================================
# Optimiser Parameters
#==============================================================================
n_outer_step = 5

# branching limit for global pruning during ply orientation optimisation
global_node_limit = 10
# branching limit for local pruning during ply orientation optimisation
local_node_limit = 10
# branching limit for global pruning at the penultimate level during ply
# orientation optimisation
global_node_limit_p = 10
# branching limit for local pruning at the last level during ply
# orientation optimisation
local_node_limit_final = 1

### Techniques to enforce the constraints
# repair to improve the convergence towards the in-plane lamination parameter
# targets
repair_membrane_switch = True
# repair to improve the convergence towards the out-of-plane lamination
# parameter targets
repair_flexural_switch = True

# penalty for the 10% rule based on ply count restrictions
penalty_10_pc_switch = False
# penalty for the 10% rule based on lamination parameter restrictions
penalty_10_lampam_switch = False
# penalty for in-plane orthotropy, based on lamination parameters
penalty_ipo_switch = False
# penalty for balance, based on ply counts
penalty_bal_switch = False
# balanced laminate scheme
balanced_scheme = False

# Coefficient for the 10% rule penalty
coeff_10 = 1
# Coefficients for the in-plane orthotropy penalty or the balance penalty
coeff_bal_ipo = 1
# Coefficient for the out-of-plane orthotropy penalty
coeff_oopo = 1

# percentage of laminate thickness for plies that can be modified during
# the refinement of membrane properties
p_A = 80
# number of plies in the last permutation during repair for disorientation
# and/or contiguity
n_D1 = 6
# number of ply shifts tested at each step of the re-designing process during
# refinement of flexural properties
n_D2 = 10
# number of times the algorithms 1 and 2 are repeated during the flexural
# property refinement
n_D3 = 2

### Other parameters
optimisation_type = 'AD'

# Lamination parameters to be considered in the multi-objective functions
lampam_to_be_optimised = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

# Lamination parameters sensitivities from the first-lebel optimiser
first_level_sensitivities = np.ones((12,), float)

# Minimum group size allowed for the smallest groups
group_size_min = 5
# Desired number of plies for the groups at each outer loop
group_size_max = np.array([1000, 12, 12, 12, 12])

#==============================================================================
timesK = [3600, 7.10, 7.40, 0.95, 1.40, 5.87, 5.30, 0.64, 0.91,
                          3600, 5.95, 6.34, 0.90, 1.59, 9.06, 6.55, 1.01, 1.37
                          ]

comments = ['Normal',
            'Symmetric',
            'Sym + contiguity',
            'Sym + disorientation',
            'Sym + disorientation',
            'Sym + 10%',
            'Sym + 10% + contiguity + damtol',
            'Sym + 10% + contiguity + damtol + disorientation',
            'Sym + 10% + contiguity + damtol + disorientation',
            'Balanced',
            'Sym + bal',
            'Sym + bal + contiguity',
            'Sym + bal + disorientation',
            'Sym + bal + disorientation',
            'Sym + bal + 10%',
            'Sym + bal + 10% + contiguity + damtol',
            'Sym + bal + 10% + contiguity + damtol + disorientation',
            'Sym + bal + 10% + contiguity + damtol + disorientation',
            ]

ply_counts = [28, 28, 28, 28, 29, 28, 28, 28, 29,
              28, 28, 28, 28, 29, 28, 28, 28, 29
              ]

lampam_targets = [
    np.array([-0.168, -0.0854, 0.0097, 0,
         -0.0072, -0.0072, -0.0072, 0,
         0.0746, -0.7087, -0.0261, 0]),
    np.array([-0.1913, -0.0612, -0.0344, 0,
         0, 0, 0, 0,
         0.0259, -0.7922, -0.0303, 0]),
    np.array([-0.0888, -0.2551, 0.0856, 0,
         0, 0, 0, 0,
         0.0628, -0.8113, -0.0123, 0]),
    np.array([-0.1542, -0.0802, 0, 0,
         -0.029, -0.029, -0.029, 0,
         0.0299, -0.8037, -0.0598, 0]),
    np.array([-0.1519, -0.0621, 0, 0,
         0, 0, 0, 0,
         0.0437, -0.79, -0.0233, 0]),
    np.array([-0.1196, -0.0585, 0, 0,
         0, 0, 0, 0,
         0.0483, -0.721, -0.0196, 0])
    ]

ssK = [
        np.array([-45, -45, 45, -45, 45, 45, 0, 0, 45, 90, 90, 90, 90, 45, 90, 90, 90, 90, 90, 0, -45, 0, -45, 45, -45, -45, 45, 45], int),
        np.array([45, -45, 45, -45, -45, 45, -45, 0, 0, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 0, 0, -45, 45, -45, -45, 45, -45, 45], int),
        np.array([45, -45, -45, 45, -45, 45, 0, -45, 45, 90, 90, 0, 90, 90, 90, 90, 0, 90, 90, 45, -45, 0, 45, -45, 45, -45, -45, 45 ], int),
        np.array([-45, -45, 0, 45, 45, 45, 45, 90, -45, 90, 90, 90, -45, 0, 0, -45, 90, 90, 90, -45, 90, 45, 45, 45, 45, 0, -45, -45], int),
        np.array([-45, -45, 0, 45, 45, 45, 45, 90, -45, 90, 90, 90, -45, 0, 0, 0, -45, 90, 90, 90, -45, 90,  45, 45,  45, 45, 0, -45, -45], int),

        np.array([45, -45, -45, -45, 45, 45, 45, 0, 0, -45, 90, 45, 90, 90, 90, 90, 45, 90, -45, 0, 0, 45, 45, 45, -45, -45, -45, 45], int),
        np.array([45, -45, -45, -45, 45, 45, 45, 0, 0, -45, 90, 45, 90, 90, 90, 90, 45, 90, -45, 0, 0, 45, 45, 45, -45, -45, -45, 45], int),
        np.array([-45, -45, 0, 45, 45, 45, 45, 90, 45, 90, -45, 90, -45, 0, 0, -45, 90, -45, 90, 45, 90, 45, 45, 45, 45, 0, -45, -45], int),
        np.array([-45, -45, 0, 45, 45, 45, 45, 90, -45, 90, 45, 0, -45, 90, 90, 90, -45, 0, 45, 90, -45, 90, 45, 45, 45, 45, 0, -45, -45], int),
        np.array([-45, 45, 45, -45, 45, -45, 0, -45, 0, 45, 90, 90, 90, 90, 90, 90, 90, 0, 0, 90, 45, 45, -45, 45, -45, -45, -45, 45], int),

        np.array([-45, 45, -45, 45, 45, -45, 0, 45, -45, 90, 0, 90, 90, 90, 90, 90, 90, 0, 90, -45, 45, 0, -45, 45, 45, -45, 45, -45], int),
        np.array([-45, 45, -45, 45, 45, -45, 0, 45, -45, 90, 90, 0, 90, 90, 90, 90, 0, 90, 90, -45, 45, 0, -45, 45, 45, -45, 45, -45], int),
        np.array([-45, -45, 0, 45, 45, 45, 45, 90, -45, 90, 90, 90, -45, 0, 0, -45, 90, 90, 90, -45, 90, 45, 45, 45, 45, 0, -45, -45], int),
        np.array([-45, -45, 0, 45, 45, 45, 45, 90, -45, 90, 90, 90, -45, 0, 0, 0, -45, 90, 90, 90, -45, 90, 45, 45, 45, 45, 0, -45, -45], int),
        np.array([45, -45, -45, 45, 45, -45, 0, -45, 0, 90, 90, 90, 90, 45, 45, 90, 90, 90, 90, 0, -45, 0, -45, 45, 45, -45, -45, 45], int),

        np.array([45, -45, -45, 45, 45, -45, 0, -45, 0, 90, 90, 90, 90, 45, 45, 90, 90, 90, 90, 0, -45, 0, -45, 45, 45, -45, -45, 45], int),
        np.array([-45, -45, 0, 45, 45, 45, 45, 90, -45, 90, 90, 90, -45, 0, 0, -45, 90, 90, 90, -45, 90, 45, 45, 45, 45, 0, -45, -45], int),
        np.array([-45, -45, 0, 45, 45, 45, 45, 90, -45, 90, 90, 90, -45, 0, 0, 0, -45, 90, 90, 90, -45, 90, 45, 45, 45, 45, 0, -45, -45], int),
        ]

for i in range(18):

    print('\n ipop', i)

    comment = comments[i]
    time_global_layerwise = timesK[i]
    n_plies = ply_counts[i]
    lampamK = calc_lampam(ssK[i])

    if i == 0: # 'Normal'
        lampam_target = lampam_targets[0]
        constraints = Constraints(
            sym=False,
            bal=False,
            dam_tol=False,
            diso=False,
            contig=False,
            rule_10_percent=False,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 1: # 'Symmetric'
        lampam_target = lampam_targets[1]
        constraints = Constraints(
            sym=True,
            bal=False,
            dam_tol=False,
            diso=False,
            contig=False,
            rule_10_percent=False,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 2: # 'Sym + contiguity'
        lampam_target = lampam_targets[1]
        constraints = Constraints(
            sym=True,
            bal=False,
            dam_tol=False,
            diso=False,
            contig=True,
            rule_10_percent=False,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 3 or i == 4: # 'Sym + disorientation'
        lampam_target = lampam_targets[1]
        constraints = Constraints(
            sym=True,
            bal=False,
            dam_tol=False,
            diso=True,
            contig=False,
            rule_10_percent=False,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 5: # 'Sym + 10%'
        lampam_target = lampam_targets[2]
        constraints = Constraints(
            sym=True,
            bal=False,
            dam_tol=False,
            diso=False,
            contig=False,
            rule_10_percent=True,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 6: # 'Sym + 10% + contiguity + damtol'
        lampam_target = lampam_targets[2]
        constraints = Constraints(
            sym=True,
            bal=False,
            dam_tol=True,
            diso=False,
            contig=True,
            rule_10_percent=True,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 7 or i == 8: # 'Sym + 10% + contiguity + damtol + disorientation'
        lampam_target = lampam_targets[2]
        constraints = Constraints(
            sym=True,
            bal=False,
            dam_tol=True,
            diso=True,
            contig=True,
            rule_10_percent=True,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 9: # 'Balanced'
        lampam_target = lampam_targets[3]
        constraints = Constraints(
            sym=False,
            bal=True,
            dam_tol=False,
            diso=False,
            contig=False,
            rule_10_percent=False,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 10: # 'Sym + bal'
        lampam_target = lampam_targets[4]
        constraints = Constraints(
            sym=True,
            bal=True,
            dam_tol=False,
            diso=False,
            contig=False,
            rule_10_percent=False,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 11: # 'Sym + bal + contiguity'
        lampam_target = lampam_targets[4]
        constraints = Constraints(
            sym=True,
            bal=True,
            dam_tol=False,
            diso=False,
            contig=True,
            rule_10_percent=False,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 12 or i == 13: # 'Sym + bal + disorientation'
        lampam_target = lampam_targets[4]
        constraints = Constraints(
            sym=True,
            bal=True,
            dam_tol=False,
            diso=True,
            contig=False,
            rule_10_percent=False,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 14: # 'Sym + bal + 10%'
        lampam_target = lampam_targets[5]
        constraints = Constraints(
            sym=True,
            bal=True,
            dam_tol=False,
            diso=False,
            contig=False,
            rule_10_percent=True,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 15: # 'Sym + bal + 10% + contiguity + damtol'
        lampam_target = lampam_targets[5]
        constraints = Constraints(
            sym=True,
            bal=True,
            dam_tol=True,
            diso=False,
            contig=True,
            rule_10_percent=True,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))
    elif i == 16 or i == 17: # 'Sym + bal + 10% + contiguity + damtol + diso'
        lampam_target = lampam_targets[5]
        constraints = Constraints(
            sym=True,
            bal=True,
            dam_tol=True,
            diso=True,
            contig=True,
            rule_10_percent=True,
            percent_0=percent_0,
            percent_45=percent_45,
            percent_90=percent_90,
            percent_135=percent_135,
            percent_45_135=percent_45_135,
            combine_45_135=combine_45_135,
            n_contig=n_contig,
            delta_angle=delta_angle,
            dam_tol_rule=dam_tol_rule,
            set_of_angles=np.array([-45, 0, 45, 90], int))

    targets = Targets(n_plies=ply_counts[i], lampam=lampam_target)
    parameters = Parameters(
        constraints=constraints,
        coeff_10=coeff_10,
        coeff_bal_ipo=coeff_bal_ipo,
        coeff_oopo=coeff_oopo,
        p_A=p_A,
        n_D1=n_D1,
        n_D2=n_D2,
        n_D3=n_D3,
        n_outer_step=n_outer_step,
        group_size_min=group_size_min,
        group_size_max=group_size_max,
        first_level_sensitivities=first_level_sensitivities,
        lampam_to_be_optimised=lampam_to_be_optimised,
        global_node_limit=global_node_limit,
        local_node_limit=local_node_limit,
        global_node_limit_p=global_node_limit_p,
        local_node_limit_final=local_node_limit_final,
        repair_membrane_switch=repair_membrane_switch,
        repair_flexural_switch=repair_flexural_switch,
        penalty_10_lampam_switch=penalty_10_lampam_switch,
        penalty_10_pc_switch=penalty_10_pc_switch,
        penalty_ipo_switch=penalty_ipo_switch,
        penalty_bal_switch=penalty_bal_switch,
        type_obj_func=1)
    #print(parameters)

    t = time.time()
    result = LAYLA_optimiser(parameters, constraints, targets, mat_prop)
    elapsed1 = time.time() - t


    check_lay_up_rules(result.ss, constraints)
    check_lay_up_rules(ssK[i], constraints)


    print('Time', elapsed1)
    print('objective with modified lamination parameter weightings',
          result.objective)

    table_result = pd.DataFrame()

#    # number of the outer loop with the best results
#    table_result.loc[0, 'best outer loop'] \
#    = result.n_outer_step_best_solution
#
#    # Number of iterations
#    table_result.loc[0, 'n_outer_step_performed'] \
#    = result.number_of_outer_steps_performed

    table_result.loc[0, 'constraints'] = comments[i]

    # Laminate ply count
    table_result.loc[0 + 0, 'Ply count'] = np.NaN
    table_result.loc[0 + 1, 'Ply count'] = ssK[i].size
    table_result.loc[0 + 2, 'Ply count'] = result.ss.size

    # objective
    table_result.loc[
        0 + 0, 'objective with initial lamination parameter weightings'] \
        = np.NaN
    table_result.loc[
        0 + 1, 'objective with initial lamination parameter weightings'] \
        = objectives(
            lampam=lampamK,
            targets=targets,
            lampam_weightings=parameters.lampam_weightings_ini,
            constraints=constraints,
            parameters=parameters)
    table_result.loc[
        0 + 2, 'objective with initial lamination parameter weightings'] \
        = objectives(
            lampam=result.lampam,
            targets=targets,
            lampam_weightings=parameters.lampam_weightings_ini,
            constraints=constraints,
            parameters=parameters)

    table_result.loc[
        0 + 0, 'objective with modified lamination parameter weightings'] \
        = np.NaN
    table_result.loc[
        0 + 1, 'objective with modified lamination parameter weightings'] \
        = objectives(
            lampam=lampamK,
            targets=targets,
            lampam_weightings=parameters.lampam_weightings_final,
            constraints=constraints,
            parameters=parameters)
    table_result.loc[
        0 + 2, 'objective with modified lamination parameter weightings'] \
        = result.objective

#    # Inhomogeneity factor
#    table_result.loc[0, 'target inhomogeneity factor'] \
#    = np.linalg.norm(lampam_target[0:4] - lampam_target[8:12])
#
#    # objectives
#    for k in range(parameters.n_outer_step):
#        table_result.loc[
#            0, f'objective iteration {k+1}'] = result.obj_tab[k]

    # Computational time in s
    table_result.loc[0 + 0, 'time (s)'] = np.NaN
    table_result.loc[0 + 1, 'time (s)'] = timesK[i]
    table_result.loc[0 + 2, 'time (s)'] = elapsed1

    # lampam_target
    table_result.loc[0 + 0, 'information'] = 'lampam_target'
    table_result.loc[0 + 0, 'index 1'] = lampam_target[0]
    table_result.loc[0 + 0, 'index 2'] = lampam_target[1]
    table_result.loc[0 + 0, 'index 3'] = lampam_target[2]
    table_result.loc[0 + 0, 'index 4'] = lampam_target[3]
    table_result.loc[0 + 0, 'index 5'] = lampam_target[4]
    table_result.loc[0 + 0, 'index 6'] = lampam_target[5]
    table_result.loc[0 + 0, 'index 7'] = lampam_target[6]
    table_result.loc[0 + 0, 'index 8'] = lampam_target[7]
    table_result.loc[0 + 0, 'index 9'] = lampam_target[8]
    table_result.loc[0 + 0, 'index 10'] = lampam_target[9]
    table_result.loc[0 + 0, 'index 11'] = lampam_target[10]
    table_result.loc[0 + 0, 'index 12'] = lampam_target[11]

    table_result.loc[0 + 1, 'information'] = 'LP errors Kennedy'
    table_result.loc[0 + 1, 'index 1'] \
    = abs(lampam_target[0] - lampamK[0])
    table_result.loc[0 + 1, 'index 2'] \
    = abs(lampam_target[1] - lampamK[1])
    table_result.loc[0 + 1, 'index 3'] \
    = abs(lampam_target[2] - lampamK[2])
    table_result.loc[0 + 1, 'index 4'] \
    = abs(lampam_target[3] - lampamK[3])
    table_result.loc[0 + 1, 'index 5'] \
    = abs(lampam_target[4] - lampamK[4])
    table_result.loc[0 + 1, 'index 6'] \
    = abs(lampam_target[5] - lampamK[5])
    table_result.loc[0 + 1, 'index 7'] \
    = abs(lampam_target[6] - lampamK[6])
    table_result.loc[0 + 1, 'index 8'] \
    = abs(lampam_target[7] - lampamK[7])
    table_result.loc[0 + 1, 'index 9'] \
    = abs(lampam_target[8] - lampamK[8])
    table_result.loc[0 + 1, 'index 10'] \
    = abs(lampam_target[9] - lampamK[9])
    table_result.loc[0 + 1, 'index 11'] \
    = abs(lampam_target[10] - lampamK[10])
    table_result.loc[0 + 1, 'index 12'] \
    = abs(lampam_target[11] - lampamK[11])

    table_result.loc[0 + 2, 'information'] = 'LP errors LAYLA'
    table_result.loc[0 + 2, 'index 1'] \
    = abs(lampam_target[0] - result.lampam[0])
    table_result.loc[0 + 2, 'index 2'] \
    = abs(lampam_target[1] - result.lampam[1])
    table_result.loc[0 + 2, 'index 3'] \
    = abs(lampam_target[2] - result.lampam[2])
    table_result.loc[0 + 2, 'index 4'] \
    = abs(lampam_target[3] - result.lampam[3])
    table_result.loc[0 + 2, 'index 5'] \
    = abs(lampam_target[4] - result.lampam[4])
    table_result.loc[0 + 2, 'index 6'] \
    = abs(lampam_target[5] - result.lampam[5])
    table_result.loc[0 + 2, 'index 7'] \
    = abs(lampam_target[6] - result.lampam[6])
    table_result.loc[0 + 2, 'index 8'] \
    = abs(lampam_target[7] - result.lampam[7])
    table_result.loc[0 + 2, 'index 9'] \
    = abs(lampam_target[8] - result.lampam[8])
    table_result.loc[0 + 2, 'index 10'] \
    = abs(lampam_target[9] - result.lampam[9])
    table_result.loc[0 + 2, 'index 11'] \
    = abs(lampam_target[10] - result.lampam[10])
    table_result.loc[0 + 2, 'index 12'] \
    = abs(lampam_target[11] - result.lampam[11])

#    # Retrieved stacking sequence at step 1
#    ss_flatten = np.array(result.ss_tab[0], dtype=str)
#    ss_flatten = ' '.join(ss_flatten)
#    table_result.loc[0, 'ss retrieved at step 1'] = ss_flatten

    # Retrieved stacking sequence
    table_result.loc[0 + 0, 'stacking sequence'] = np.NaN
    ss_flattenK = np.array(ssK[i], dtype=str)
    ss_flattenK = ' '.join(ss_flattenK)
    table_result.loc[0 + 1, 'stacking sequence'] = ss_flattenK
    ss_flatten = np.array(result.ss, dtype=str)
    ss_flatten = ' '.join(ss_flatten)
    table_result.loc[0 + 2, 'stacking sequence'] = ss_flatten


#    # Target stacking sequence
#    table_result.loc[0, 'ss target'] = np.NaN

#    # Ply counts
#    table_result.loc[0, 'N0_target'] = np.NaN
#    table_result.loc[0, 'N90_target'] = np.NaN
#    table_result.loc[0, 'N45_target'] = np.NaN
#    table_result.loc[0, 'N-45_target'] = np.NaN
#    N0 = sum(result.ss == 0)
#    N90 = sum(result.ss == 90)
#    N45 = sum(result.ss == 45)
#    N135 = sum(result.ss == -45)
#    table_result.loc[0, 'N0 - N0_target'] = np.NaN
#    table_result.loc[0, 'N90 - N90_target'] = np.NaN
#    table_result.loc[0, 'N45 - N45_target'] = np.NaN
#    table_result.loc[0, 'N-45 - N-45_target'] = np.NaN
#    table_result.loc[0, 'penalty value for the 10% rule'] \
#    = calc_penalty_10_ss(result.ss, constraints)

#    for ind in range(n_outer_step):
#        # numbers of stacks at the last level of the last group search
#        table_result.loc[0, 'n_designs_last_level ' + str(ind + 1)] \
#        = result.n_designs_last_level_tab[ind]
#        # numbers of repaired stacks at the last group search
#        table_result.loc[0, 'n_designs_repaired ' + str(ind + 1)] \
#        = result.n_designs_repaired_tab[ind]
#        # numbers of unique repaired stacks at the last group search
#        table_result.loc[0, 'n_designs_repaired_unique ' + str(ind + 1)] \
#        = result.n_designs_repaired_unique_tab[ind]

#    # in-plane orthotropy
#    ipo_now = ipo_param_1_12(result.lampam, mat_prop, constraints.sym)
#    table_result.loc[0, 'In-plane orthotropy parameter 1'] = ipo_now[0]
#    table_result.loc[0, 'In-plane orthotropy parameter 2'] = ipo_now[1]
#    table_result.loc[0, 'In-plane orthotropy parameter 3'] = ipo_now[2]
#    table_result.loc[0, 'In-plane orthotropy parameter 4'] = ipo_now[3]
#    table_result.loc[0, 'In-plane orthotropy parameter 5'] = ipo_now[4]
#    table_result.loc[0, 'In-plane orthotropy parameter 6'] = ipo_now[5]
#    table_result.loc[0, 'In-plane orthotropy parameter 7'] = ipo_now[6]
#    table_result.loc[0, 'In-plane orthotropy parameter 8'] = ipo_now[7]
#    table_result.loc[0, 'In-plane orthotropy parameter 9'] = ipo_now[8]
#    table_result.loc[0, 'In-plane orthotropy parameter 10'] = ipo_now[9]
#    table_result.loc[0, 'In-plane orthotropy parameter 11'] = ipo_now[10]
#    table_result.loc[0, 'In-plane orthotropy parameter 12'] = ipo_now[11]

    AK = A_from_lampam(lampamK, mat_prop)
    A11K = AK[0, 0]
    A22K = AK[1, 1]
    A12K = AK[0, 1]
    A66K = AK[2, 2]
    A16K = AK[0, 2]
    A26K = AK[1, 2]

    BK = B_from_lampam(lampamK, mat_prop)
    B11K = BK[0, 0]
    B22K = BK[1, 1]
    B12K = BK[0, 1]
    B66K = BK[2, 2]
    B16K = BK[0, 2]
    B26K = BK[1, 2]

    DK = D_from_lampam(lampamK, mat_prop)
    D11K = DK[0, 0]
    D22K = DK[1, 1]
    D12K = DK[0, 1]
    D66K = DK[2, 2]
    D16K = DK[0, 2]
    D26K = DK[1, 2]

    A = A_from_lampam(result.lampam, mat_prop)
    A11 = A[0, 0]
    A22 = A[1, 1]
    A12 = A[0, 1]
    A66 = A[2, 2]
    A16 = A[0, 2]
    A26 = A[1, 2]

    B = B_from_lampam(result.lampam, mat_prop)
    B11 = B[0, 0]
    B22 = B[1, 1]
    B12 = B[0, 1]
    B66 = B[2, 2]
    B16 = B[0, 2]
    B26 = B[1, 2]

    D = D_from_lampam(result.lampam, mat_prop)
    D11 = D[0, 0]
    D22 = D[1, 1]
    D12 = D[0, 1]
    D66 = D[2, 2]
    D16 = D[0, 2]
    D26 = D[1, 2]

    A_target = A_from_lampam(lampam_target, mat_prop)
    A11_target = A_target[0, 0]
    A22_target = A_target[1, 1]
    A12_target = A_target[0, 1]
    A66_target = A_target[2, 2]
    A16_target = A_target[0, 2]
    A26_target = A_target[1, 2]

    B_target = B_from_lampam(lampam_target, mat_prop)
    B11_target = B_target[0, 0]
    B22_target = B_target[1, 1]
    B12_target = B_target[0, 1]
    B66_target = B_target[2, 2]
    B16_target = B_target[0, 2]
    B26_target = B_target[1, 2]

    D_target = D_from_lampam(lampam_target, mat_prop)
    D11_target = D_target[0, 0]
    D22_target = D_target[1, 1]
    D12_target = D_target[0, 1]
    D66_target = D_target[2, 2]
    D16_target = D_target[0, 2]
    D26_target = D_target[1, 2]

    table_result.loc[0 + 0, 'stifnesses'] = 'targets'
    table_result.loc[0 + 0, 'A11'] = A11_target
    table_result.loc[0 + 0, 'A22'] = A22_target
    table_result.loc[0 + 0, 'A12'] = A12_target
    table_result.loc[0 + 0, 'A66'] = A66_target
    table_result.loc[0 + 0, 'A16'] = A16_target
    table_result.loc[0 + 0, 'A26'] = A26_target

    table_result.loc[0 + 0, 'B11'] = B11_target
    table_result.loc[0 + 0, 'B22'] = B22_target
    table_result.loc[0 + 0, 'B12'] = B12_target
    table_result.loc[0 + 0, 'B66'] = B66_target
    table_result.loc[0 + 0, 'B16'] = B16_target
    table_result.loc[0 + 0, 'B26'] = B26_target

    table_result.loc[0 + 0, 'D11'] = D11_target
    table_result.loc[0 + 0, 'D22'] = D22_target
    table_result.loc[0 + 0, 'D12'] = D12_target
    table_result.loc[0 + 0, 'D66'] = D66_target
    table_result.loc[0 + 0, 'D16'] = D16_target
    table_result.loc[0 + 0, 'D26'] = D26_target

    table_result.loc[0 + 1, 'stifnesses'] = 'Kennedy'
    table_result.loc[0 + 1, 'A11'] = A11K
    table_result.loc[0 + 1, 'A22'] = A22K
    table_result.loc[0 + 1, 'A12'] = A12K
    table_result.loc[0 + 1, 'A66'] = A66K
    table_result.loc[0 + 1, 'A16'] = A16K
    table_result.loc[0 + 1, 'A26'] = A26K

    table_result.loc[0 + 1, 'B11'] = B11K
    table_result.loc[0 + 1, 'B22'] = B22K
    table_result.loc[0 + 1, 'B12'] = B12K
    table_result.loc[0 + 1, 'B66'] = B66K
    table_result.loc[0 + 1, 'B16'] = B16K
    table_result.loc[0 + 1, 'B26'] = B26K

    table_result.loc[0 + 1, 'D11'] = D11K
    table_result.loc[0 + 1, 'D22'] = D22K
    table_result.loc[0 + 1, 'D12'] = D12K
    table_result.loc[0 + 1, 'D66'] = D66K
    table_result.loc[0 + 1, 'D16'] = D16K
    table_result.loc[0 + 1, 'D26'] = D26K

    table_result.loc[0 + 2, 'stifnesses'] = 'LAYLA'
    table_result.loc[0 + 2, 'A11'] = A11
    table_result.loc[0 + 2, 'A22'] = A22
    table_result.loc[0 + 2, 'A12'] = A12
    table_result.loc[0 + 2, 'A66'] = A66
    table_result.loc[0 + 2, 'A16'] = A16
    table_result.loc[0 + 2, 'A26'] = A26

    table_result.loc[0 + 2, 'B11'] = B11
    table_result.loc[0 + 2, 'B22'] = B22
    table_result.loc[0 + 2, 'B12'] = B12
    table_result.loc[0 + 2, 'B66'] = B66
    table_result.loc[0 + 2, 'B16'] = B16
    table_result.loc[0 + 2, 'B26'] = B26

    table_result.loc[0 + 2, 'D11'] = D11
    table_result.loc[0 + 2, 'D22'] = D22
    table_result.loc[0 + 2, 'D12'] = D12
    table_result.loc[0 + 2, 'D66'] = D66
    table_result.loc[0 + 2, 'D16'] = D16
    table_result.loc[0 + 2, 'D26'] = D26

#    if A11_target:
#        table_result.loc[0 + 1, 'diff A11 percentage']\
#        = 100 * abs((A11K - A11_target)/A11_target)
#    else:
#        table_result.loc[0 + 1, 'diff A11 percentage'] = 0
#    if A22_target:
#        table_result.loc[0 + 1, 'diff A22 percentage']\
#        = 100 * abs((A22K - A22_target)/A22_target)
#    else:
#        table_result.loc[0 + 1, 'diff A22 percentage'] = 0
#    if A12_target:
#        table_result.loc[0 + 1, 'diff A12 percentage']\
#        = 100 * abs((A12K - A12_target)/A12_target)
#    else:
#        table_result.loc[0 + 1, 'diff A12 percentage'] = 0
#    if A66_target:
#        table_result.loc[0 + 1, 'diff A66 percentage']\
#        = 100 * abs((A66K - A66_target)/A66_target)
#    else:
#        table_result.loc[0 + 1, 'diff A66 percentage'] = 0
#    if A16_target:
#        table_result.loc[0 + 1, 'diff A16 percentage']\
#        = 100 * abs((A16K - A16_target)/A16_target)
#    else:
#        table_result.loc[0 + 1, 'diff A16 percentage'] = 0
#    if A26_target:
#        table_result.loc[0 + 1, 'diff A26 percentage']\
#        = 100 * abs((A26K - A26_target)/A26_target)
#    else:
#        table_result.loc[0 + 1, 'diff A26 percentage'] = 0
#
#    if B11_target:
#        table_result.loc[0 + 1, 'diff B11 percentage']\
#        = 100 * abs((B11K - B11_target)/B11_target)
#    else:
#        table_result.loc[0 + 1, 'diff B11 percentage'] = 0
#    if B22_target:
#        table_result.loc[0 + 1, 'diff B22 percentage']\
#        = 100 * abs((B22K - B22_target)/B22_target)
#    else:
#        table_result.loc[0 + 1, 'diff B22 percentage'] = 0
#    if B12_target:
#        table_result.loc[0 + 1, 'diff B12 percentage']\
#        = 100 * abs((B12K - B12_target)/B12_target)
#    else:
#        table_result.loc[0 + 1, 'diff B12 percentage'] = 0
#    if B66_target:
#        table_result.loc[0 + 1, 'diff B66 percentage']\
#        = 100 * abs((B66K - B66_target)/B66_target)
#    else:
#        table_result.loc[0 + 1, 'diff B66 percentage'] = 0
#    if B16_target:
#        table_result.loc[0 + 1, 'diff B16 percentage']\
#        = 100 * abs((B16K - B16_target)/B16_target)
#    else:
#        table_result.loc[0 + 1, 'diff B16 percentage'] = 0
#    if B26_target:
#        table_result.loc[0 + 1, 'diff B26 percentage']\
#        = 100 * abs((B26K - B26_target)/B26_target)
#    else:
#        table_result.loc[0 + 1, 'diff B26 percentage'] = 0
#
#    if D11_target:
#        table_result.loc[0 + 1, 'diff D11 percentage']\
#        = 100 * abs((D11K - D11_target)/D11_target)
#    else:
#        table_result.loc[0 + 1, 'diff D11 percentage'] = 0
#    if D22_target:
#        table_result.loc[0 + 1, 'diff D22 percentage']\
#        = 100 * abs((D22K - D22_target)/D22_target)
#    else:
#        table_result.loc[0 + 1, 'diff D22 percentage'] = 0
#    if D12_target:
#        table_result.loc[0 + 1, 'diff D12 percentage']\
#        = 100 * abs((D12K - D12_target)/D12_target)
#    else:
#        table_result.loc[0 + 1, 'diff D12 percentage'] = 0
#    if D66_target:
#        table_result.loc[0 + 1, 'diff D66 percentage']\
#        = 100 * abs((D66K - D66_target)/D66_target)
#    else:
#        table_result.loc[0 + 1, 'diff D66 percentage'] = 0
#    if D16_target:
#        table_result.loc[0 + 1, 'diff D16 percentage']\
#        = 100 * abs((D16K - D16_target)/D16_target)
#    else:
#        table_result.loc[0 + 1, 'diff D16 percentage'] = 0
#    if D26_target:
#        table_result.loc[0 + 1, 'diff D26 percentage']\
#        = 100 * abs((D26K - D26_target)/D26_target)
#    else:
#        table_result.loc[0 + 1, 'diff D26 percentage'] = 0
#
#
#
#
#    if A11_target:
#        table_result.loc[0 + 2, 'diff A11 percentage']\
#        = 100 * abs((A11 - A11_target)/A11_target)
#    else:
#        table_result.loc[0 + 2, 'diff A11 percentage'] = 0
#    if A22_target:
#        table_result.loc[0 + 2, 'diff A22 percentage']\
#        = 100 * abs((A22 - A22_target)/A22_target)
#    else:
#        table_result.loc[0 + 2, 'diff A22 percentage'] = 0
#    if A12_target:
#        table_result.loc[0 + 2, 'diff A12 percentage']\
#        = 100 * abs((A12 - A12_target)/A12_target)
#    else:
#        table_result.loc[0 + 2, 'diff A12 percentage'] = 0
#    if A66_target:
#        table_result.loc[0 + 2, 'diff A66 percentage']\
#        = 100 * abs((A66 - A66_target)/A66_target)
#    else:
#        table_result.loc[0 + 2, 'diff A66 percentage'] = 0
#    if A16_target:
#        table_result.loc[0 + 2, 'diff A16 percentage']\
#        = 100 * abs((A16 - A16_target)/A16_target)
#    else:
#        table_result.loc[0 + 2, 'diff A16 percentage'] = 0
#    if A26_target:
#        table_result.loc[0 + 2, 'diff A26 percentage']\
#        = 100 * abs((A26 - A26_target)/A26_target)
#    else:
#        table_result.loc[0 + 2, 'diff A26 percentage'] = 0
#
#    if B11_target:
#        table_result.loc[0 + 2, 'diff B11 percentage']\
#        = 100 * abs((B11 - B11_target)/B11_target)
#    else:
#        table_result.loc[0 + 2, 'diff B11 percentage'] = 0
#    if B22_target:
#        table_result.loc[0 + 2, 'diff B22 percentage']\
#        = 100 * abs((B22 - B22_target)/B22_target)
#    else:
#        table_result.loc[0 + 2, 'diff B22 percentage'] = 0
#    if B12_target:
#        table_result.loc[0 + 2, 'diff B12 percentage']\
#        = 100 * abs((B12 - B12_target)/B12_target)
#    else:
#        table_result.loc[0 + 2, 'diff B12 percentage'] = 0
#    if B66_target:
#        table_result.loc[0 + 2, 'diff B66 percentage']\
#        = 100 * abs((B66 - B66_target)/B66_target)
#    else:
#        table_result.loc[0 + 2, 'diff B66 percentage'] = 0
#    if B16_target:
#        table_result.loc[0 + 2, 'diff B16 percentage']\
#        = 100 * abs((B16 - B16_target)/B16_target)
#    else:
#        table_result.loc[0 + 2, 'diff B16 percentage'] = 0
#    if B26_target:
#        table_result.loc[0 + 2, 'diff B26 percentage']\
#        = 100 * abs((B26 - B26_target)/B26_target)
#    else:
#        table_result.loc[0 + 2, 'diff B26 percentage'] = 0
#
#    if D11_target:
#        table_result.loc[0 + 2, 'diff D11 percentage']\
#        = 100 * abs((D11 - D11_target)/D11_target)
#    else:
#        table_result.loc[0 + 2, 'diff D11 percentage'] = 0
#    if D22_target:
#        table_result.loc[0 + 2, 'diff D22 percentage']\
#        = 100 * abs((D22 - D22_target)/D22_target)
#    else:
#        table_result.loc[0 + 2, 'diff D22 percentage'] = 0
#    if D12_target:
#        table_result.loc[0 + 2, 'diff D12 percentage']\
#        = 100 * abs((D12 - D12_target)/D12_target)
#    else:
#        table_result.loc[0 + 2, 'diff D12 percentage'] = 0
#    if D66_target:
#        table_result.loc[0 + 2, 'diff D66 percentage']\
#        = 100 * abs((D66 - D66_target)/D66_target)
#    else:
#        table_result.loc[0 + 2, 'diff D66 percentage'] = 0
#    if D16_target:
#        table_result.loc[0 + 2, 'diff D16 percentage']\
#        = 100 * abs((D16 - D16_target)/D16_target)
#    else:
#        table_result.loc[0 + 2, 'diff D16 percentage'] = 0
#    if D26_target:
#        table_result.loc[0 + 2, 'diff D26 percentage']\
#        = 100 * abs((D26 - D26_target)/D26_target)
#    else:
#        table_result.loc[0 + 2, 'diff D26 percentage'] = 0

#    print(table_result)

    append_df_to_excel(
        result_filename, table_result, 'results', index=True, header=True)


### Write results in a excell sheet
save_constraints_LAYLA(result_filename, constraints)
save_parameters_LAYLA_V02(result_filename, parameters)
save_materials(result_filename, mat_prop)
autofit_column_widths(result_filename)


