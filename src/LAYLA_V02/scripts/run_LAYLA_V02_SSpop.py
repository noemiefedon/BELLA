# -*- coding: utf-8 -*-
"""
LAYLA retrieves the LAminate LAY-ups from lamination parameters

LAYLA_sspop is a script applying the optimiser LAYLA to sets of target
lamination parameters.

These lamination parameters come from input files and are associated to the
poulations of stacking sequences created in the folder Populations.
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np
import pandas as pd
import time
import sys
sys.path.append(r'C:\BELLA_and_LAYLA')
from src.LAYLA_V02.parameters import Parameters
from src.LAYLA_V02.constraints import Constraints
from src.LAYLA_V02.targets import Targets
from src.LAYLA_V02.optimiser import LAYLA_optimiser
from src.LAYLA_V02.materials import Material
from src.LAYLA_V02.objectives import objectives
from src.guidelines.ipo_oopo import ipo_param_1_12
from src.guidelines.ipo_oopo import calc_penalty_ipo_param
from src.guidelines.ten_percent_rule import calc_penalty_10_ss
from src.CLA.ABD import A_from_lampam, B_from_lampam, D_from_lampam

from src.divers.excel import autofit_column_widths
from src.divers.excel import delete_file
from src.LAYLA_V02.save_set_up import save_constraints_LAYLA
from src.LAYLA_V02.save_set_up import save_parameters_LAYLA_V02
from src.LAYLA_V02.save_set_up import save_materials

from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

#==============================================================================
# Target population
#==============================================================================
# Number of plies
n_plies = 40
#n_plies = 80
#n_plies = 200
#==============================================================================
# Results saving
#==============================================================================
filename_end = ''
#==============================================================================
# Type of optimisations
#==============================================================================
optimisation_type = 'A'
optimisation_type = 'D'
optimisation_type = 'AD'
#==============================================================================
# design and manufacturing constraints
#==============================================================================
### Set of design and manufacturing constraints:
constraints_set = 'C0'
constraints_set = 'C1'
#     C0:   - No design and manufacturing constraints other than symmetry
#     C1:   - in-plane orthotropy enforced with penalties and repair
#           - 10% rule enforced with repair
#                   - 10% 0deg plies
#                   - 10% 90 deg plies
#                   - 5% 45deg plies
#                   - 5% -45 deg plies
#           - disorientation rule with Delta(theta) = 45 deg
#           - contiguity rule with n_contig = 5

# set of admissible fibre orientations
set_of_angles = np.array([-45, 0, 45, 90], dtype=int)
#set_of_angles = np.array([-45, 0, 45, 90, +30, -30, +60, -60], dtype=int)

# symmetry
sym = True

# balance and in-plane orthotropy requirements
if constraints_set == 'C0':
    bal = False
    ipo = False
else:
    bal = True
    ipo = True

# out-of-plane orthotropy requirements
oopo = False

# damage tolerance
# rule 1: one outer ply at + or -45 deg at laminate surfaces
# rule 2: [+45, -45] or [-45, +45] plies at laminate surfaces
# rule 3: [+45, -45], [+45, +45], [-45, -45] or [-45, +45] plies at laminate
# surfaces
dam_tol = False
dam_tol_rule = 0
#if constraints_set == 'C0':
#    dam_tol = False
#    dam_tol_rule = 0
#else:
#    dam_tol = True
#    dam_tol_rule = 1
#    dam_tol_rule = 2
##    dam_tol_rule = 3

# 10% rule
if constraints_set == 'C0':
    rule_10_percent = False
else:
    rule_10_percent = True
combine_45_135 = True
percent_0 = 10 # percentage used in the 10% rule for 0 deg plies
percent_45 =  0 # percentage used in the 10% rule for +45 deg plies
percent_90 = 10 # percentage used in the 10% rule for 90 deg plies
percent_135 = 0 # percentage used in the 10% rule for -45 deg plies
percent_45_135 = 10 # percentage used in the 10% rule for +-45 deg plies

# disorientation
if constraints_set == 'C0':
    diso = False
else:
    diso = True

# Upper bound of the variation of fibre orientation between two
# contiguous plies if the disorientation constraint is active
delta_angle = 45

# contiguity
if constraints_set == 'C0':
    contig = False
else:
    contig = True

n_contig = 5
# No more that constraints.n_contig plies with same fibre orientation should be
# next to each other if the contiguity constraint is active. The value taken
# can only be 2, 3, 4 or 5, otherwise test functions should be modified

constraints = Constraints(
    sym=sym,
    bal=bal,
    ipo=ipo,
    oopo=oopo,
    dam_tol=dam_tol,
    dam_tol_rule=dam_tol_rule,
    rule_10_percent=rule_10_percent,
    percent_0=percent_0,
    percent_45=percent_45,
    percent_90=percent_90,
    percent_135=percent_135,
    percent_45_135=percent_45_135,
    combine_45_135=combine_45_135,
    diso=diso,
    contig=contig,
    n_contig=n_contig,
    delta_angle=delta_angle,
    set_of_angles=set_of_angles)

#==============================================================================
# Material properties
#==============================================================================
# Elastic modulus in the fibre direction (Pa)
E11 = 130e9
# Elastic modulus in the transverse direction (Pa)
E22 = 9e9
# Poisson's ratio relating transverse deformation and axial loading (-)
nu12 = 0.3
# In-plane shear modulus (Pa)
G12 = 4e9
mat_prop = Material(E11 = E11, E22 = E22, G12 = G12, nu12 = nu12)

#==============================================================================
# Optimiser Parameters
#==============================================================================
# number of outer loops
n_outer_step = 5

# branching limit for global pruning during ply orientation optimisation
global_node_limit = 50
# branching limit for local pruning during ply orientation optimisation
local_node_limit = 100
# branching limit for global pruning at the penultimate level during ply
# orientation optimisation
global_node_limit_p = 50
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

if constraints_set == 'C0':
    # penalty for the 10% rule based on ply count restrictions
    penalty_10_pc_switch = False
    # penalty for the 10% rule based on lamination parameter restrictions
    penalty_10_lampam_switch = False
    # penalty for in-plane orthotropy, based on lamination parameters
    penalty_ipo_switch = False
    # penalty for balance, based on ply counts
    penalty_bal_switch = False

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

# Minimum group size allowed for the smallest groups
group_size_min = 5
# Desired number of plies for the groups at each outer loop
group_size_max = np.array([1000, 8, 8, 8, 8])

# Lamination parameters to be considered in the multi-objective functions
if optimisation_type == 'A':
    if constraints.set_of_angles is np.array([-45, 0, 45, 90], int):
        lampam_to_be_optimised = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    else:
        lampam_to_be_optimised = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
if optimisation_type == 'D':
    if constraints.set_of_angles is np.array([-45, 0, 45, 90], int):
        lampam_to_be_optimised = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0])
    else:
        lampam_to_be_optimised = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1])
if optimisation_type == 'AD':
    if constraints.set_of_angles is np.array([-45, 0, 45, 90], int):
        lampam_to_be_optimised = np.array([1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0])
    else:
        lampam_to_be_optimised = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])

# Lamination parameters sensitivities from the first-lebel optimiser
first_level_sensitivities = np.ones((12,), float)

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
    penalty_bal_switch=penalty_bal_switch)

#==============================================================================
# DO NOT CHANGE FROM THIS POINT
#==============================================================================
result_filename = constraints_set + '-' + str(n_plies) + 'plies-' \
+ optimisation_type +  filename_end + '.xlsx'
delete_file(result_filename)

### Import the target lamination parameters
if constraints_set == 'C0':
    data_filename = '/LAYLA_and_BELLA/populations/pop_sym_C0_' \
    + str(n_plies) + 'plies.xlsx'
else:
    data_filename = '/LAYLA_and_BELLA/populations/pop_sym_C1_' \
    + str(n_plies) + 'plies.xlsx'

### Import the target lamination parameters
data = pd.read_excel(data_filename, sheet_name='stacks')
if data.size == 0:
    raise Exception(
            'Oops, no population of target lamination parameters found')

### Initialisation of the result columns
table_result = pd.DataFrame()

#==========================================================================
#     Optimiser Runs
#==========================================================================
for i in range(len(data.index)):
#for i in range(0, 1):

    print('\n ipop', i)
    ### Store targets
    n_plies_lam = data.loc[i, 'ply_counts']
    # Stacking sequence considered for the 'layerwise_ss' approach
    ss_ini = 0*np.ones((n_plies_lam,), dtype=int)
    lampam_target = np.empty((12,), float)
    lampam_target[0] = data.loc[i, 'lampam[1]']
    lampam_target[1] = data.loc[i, 'lampam[2]']
    lampam_target[2] = data.loc[i, 'lampam[3]']
    lampam_target[3] = data.loc[i, 'lampam[4]']
    lampam_target[4] = data.loc[i, 'lampam[5]']
    lampam_target[5] = data.loc[i, 'lampam[6]']
    lampam_target[6] = data.loc[i, 'lampam[7]']
    lampam_target[7] = data.loc[i, 'lampam[8]']
    lampam_target[8] = data.loc[i, 'lampam[9]']
    lampam_target[9] = data.loc[i, 'lampam[10]']
    lampam_target[10] = data.loc[i, 'lampam[11]']
    lampam_target[11] = data.loc[i, 'lampam[12]']

    N0_Target = data.loc[i, 'N0']
    N90_Target = data.loc[i, 'N90']
    N45_Target = data.loc[i, 'N45']
    N135_Target = data.loc[i, 'N-45']

    ss_target = data.loc[i, 'ss']

    A11_target = data.loc[i, 'A11']
    A22_target = data.loc[i, 'A22']
    A12_target = data.loc[i, 'A12']
    A66_target = data.loc[i, 'A66']
    A16_target = data.loc[i, 'A16']
    A26_target = data.loc[i, 'A26']

    B11_target = data.loc[i, 'B11']
    B22_target = data.loc[i, 'B22']
    B12_target = data.loc[i, 'B12']
    B66_target = data.loc[i, 'B66']
    B16_target = data.loc[i, 'B16']
    B26_target = data.loc[i, 'B26']

    D11_target = data.loc[i, 'D11']
    D22_target = data.loc[i, 'D22']
    D12_target = data.loc[i, 'D12']
    D66_target = data.loc[i, 'D66']
    D16_target = data.loc[i, 'D16']
    D26_target = data.loc[i, 'D26']

    targets = Targets(
        n_plies=n_plies_lam, lampam=lampam_target, stack=ss_target)
#    print('target', ss_target)

    ### Algorithm run
    print(f'Algorithm running.')
    print('Laminate type: ', constraints.laminate_scheme)
    print('Laminate type of the target stacking sequences: ',
          constraints.laminate_scheme_test)

    t = time.time()
    result = LAYLA_optimiser(parameters, constraints, targets, mat_prop)
    elapsed1 = time.time() - t

    ### Results processing and display
    if not result.completed:

        # Laminate ply count
        table_result.loc[i, 'ply count'] = n_plies_lam

        # number of the outer loop with the best results
        table_result.loc[i, 'best outer loop'] = np.NaN

        # Computational time in s
        table_result.loc[i, 'time (s)'] = np.NaN

#        # Number of objective function evaluations
#        table_result.loc[i, 'Number of objective function evaluations'] = \
#        np.NaN

        # Number of iterations
        table_result.loc[i, 'n_outer_step_performed'] = np.NaN

        # objective
        table_result.loc[
            i, 'objective with initial lamination parameter weightings'] \
            = np.NaN
        table_result.loc[
            i, 'objective with modified lamination parameter weightings'] \
            = np.NaN

        # Inhomogeneity factor
        table_result.loc[i, 'target inhomogeneity factor'] = \
        np.linalg.norm(lampam_target[0:4] - lampam_target[8:12])

        # objectives
        for k in range(parameters.n_outer_step):
            table_result.loc[i, f'objective iteration {k+1}'] = np.NaN

        # lampam_target - lampamRetrieved
        table_result.loc[i, 'error1 = abs(lampam_target[1]-lampam[1])'] \
        = np.NaN
        table_result.loc[i, 'error2'] = np.NaN
        table_result.loc[i, 'error3'] = np.NaN
        table_result.loc[i, 'error4'] = np.NaN
        table_result.loc[i, 'error5'] = np.NaN
        table_result.loc[i, 'error6'] = np.NaN
        table_result.loc[i, 'error7'] = np.NaN
        table_result.loc[i, 'error8'] = np.NaN
        table_result.loc[i, 'error9'] = np.NaN
        table_result.loc[i, 'error10'] = np.NaN
        table_result.loc[i, 'error11'] = np.NaN
        table_result.loc[i, 'error12'] = np.NaN

        # lampam_target
        table_result.loc[i, 'lampam_target[1]'] = lampam_target[0]
        table_result.loc[i, 'lampam_target[2]'] = lampam_target[1]
        table_result.loc[i, 'lampam_target[3]'] = lampam_target[2]
        table_result.loc[i, 'lampam_target[4]'] = lampam_target[3]
        table_result.loc[i, 'lampam_target[5]'] = lampam_target[4]
        table_result.loc[i, 'lampam_target[6]'] = lampam_target[5]
        table_result.loc[i, 'lampam_target[7]'] = lampam_target[6]
        table_result.loc[i, 'lampam_target[8]'] = lampam_target[7]
        table_result.loc[i, 'lampam_target[9]'] = lampam_target[8]
        table_result.loc[i, 'lampam_target[10]'] = lampam_target[9]
        table_result.loc[i, 'lampam_target[11]'] = lampam_target[10]
        table_result.loc[i, 'lampam_target[12]'] = lampam_target[11]

        # Retrieved stacking sequence at step 1
        table_result.loc[i, 'ss retrieved at step 1'] = np.NaN

        # Retrieved stacking sequence
        table_result.loc[i, 'ss retrieved'] = np.NaN

        # Target stacking sequence
        ss_flatten = np.array(ss_target, dtype=str)
        #ss_flatten = ' '.join(ss_flatten)
        table_result.loc[i, 'ss target'] = ss_flatten

        # Ply counts
        table_result.loc[i, 'N0_target'] = N0_Target
        table_result.loc[i, 'N90_target'] = N90_Target
        table_result.loc[i, 'N45_target'] = N45_Target
        table_result.loc[i, 'N-45_target'] = N135_Target
        table_result.loc[i, 'N0 - N0_target'] = np.NaN
        table_result.loc[i, 'N90 - N90_target'] = np.NaN
        table_result.loc[i, 'N45 - N45_target'] = np.NaN
        table_result.loc[i, 'N-45 - N-45_target'] = np.NaN
        table_result.loc[i, 'penalty value for the 10% rule'] = np.NaN

        for ind in range(n_outer_step):
            # numbers of stacks at the last level of the last group search
            table_result.loc[i, 'n_designs_last_level ' + str(ind + 1)] \
            = np.NaN
            # numbers of repaired stacks at the last group search
            table_result.loc[i, 'n_designs_repaired ' + str(ind + 1)] \
            = np.NaN
            # numbers of unique repaired stacks at the last group search
            table_result.loc[i, 'n_designs_repaired_unique ' + str(ind + 1)] \
            = np.NaN

        # in-plane orthotropy
        table_result.loc[i, 'In-plane orthotropy parameter 1'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 2'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 3'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 4'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 5'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 6'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 7'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 8'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 9'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 10'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 11'] = np.NaN
        table_result.loc[i, 'In-plane orthotropy parameter 12'] = np.NaN

        table_result.loc[i, 'diff A11 percentage'] = np.NaN
        table_result.loc[i, 'diff A22 percentage'] = np.NaN
        table_result.loc[i, 'diff A12 percentage'] = np.NaN
        table_result.loc[i, 'diff A66 percentage'] = np.NaN
        table_result.loc[i, 'diff A16 percentage'] = np.NaN
        table_result.loc[i, 'diff A26 percentage'] = np.NaN

        table_result.loc[i, 'diff B11 percentage'] = np.NaN
        table_result.loc[i, 'diff B22 percentage'] = np.NaN
        table_result.loc[i, 'diff B12 percentage'] = np.NaN
        table_result.loc[i, 'diff B66 percentage'] = np.NaN
        table_result.loc[i, 'diff B16 percentage'] = np.NaN
        table_result.loc[i, 'diff B26 percentage'] = np.NaN

        table_result.loc[i, 'diff D11 percentage'] = np.NaN
        table_result.loc[i, 'diff D22 percentage'] = np.NaN
        table_result.loc[i, 'diff D12 percentage'] = np.NaN
        table_result.loc[i, 'diff D66 percentage'] = np.NaN
        table_result.loc[i, 'diff D16 percentage'] = np.NaN
        table_result.loc[i, 'diff D26 percentage'] = np.NaN

#        table_result.loc[i, 'diff A11 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff A22 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff A12 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff A66 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff A16 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff A26 percentage - approx'] = np.NaN
#
#        table_result.loc[i, 'diff D11 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff D22 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff D12 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff D66 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff D16 percentage - approx'] = np.NaN
#        table_result.loc[i, 'diff D26 percentage - approx'] = np.NaN
    else:

        print('Time', elapsed1)
        print('objective with modified lamination parameter weightings',
              result.objective)

        # Laminate ply count
        table_result.loc[i, 'Ply count'] = n_plies_lam

        # number of the outer loop with the best results
        table_result.loc[i, 'best outer loop'] \
        = result.n_outer_step_best_solution

        # Computational time in s
        table_result.loc[i, 'time (s)'] = elapsed1

#        # Number of objective function evaluations
#        table_result.loc[i, 'Number of objective function evaluations'] \
#        = " ".join(result.n_obj_func_calls_tab.astype(str))

        # Number of iterations
        table_result.loc[i, 'n_outer_step_performed'] \
        = result.number_of_outer_steps_performed

        # objective
        table_result.loc[
            i, 'objective with initial lamination parameter weightings'] \
            = objectives(
                lampam=result.lampam,
                targets=targets,
                lampam_weightings=parameters.lampam_weightings_ini,
                constraints=constraints,
                parameters=parameters)

        table_result.loc[
            i, 'objective with modified lamination parameter weightings'] \
            = result.objective

        # Inhomogeneity factor
        table_result.loc[i, 'target inhomogeneity factor'] \
        = np.linalg.norm(lampam_target[0:4] - lampam_target[8:12])

        # objectives
        for k in range(parameters.n_outer_step):
            table_result.loc[
                i, f'objective iteration {k+1}'] = result.obj_tab[k]

        # lampam_target - lampamRetrieved
        table_result.loc[i, 'error1 = abs(lampam_target[1]-lampam[1])'] \
        = abs(lampam_target[0] - result.lampam[0])
        table_result.loc[i, 'error2'] = abs(
            lampam_target[1] - result.lampam[1])
        table_result.loc[i, 'error3'] = abs(
            lampam_target[2]- result.lampam[2])
        table_result.loc[i, 'error4'] = abs(
            lampam_target[3]- result.lampam[3])
        table_result.loc[i, 'error5'] = abs(
            lampam_target[4]- result.lampam[4])
        table_result.loc[i, 'error6'] = abs(
            lampam_target[5]- result.lampam[5])
        table_result.loc[i, 'error7'] = abs(
            lampam_target[6]- result.lampam[6])
        table_result.loc[i, 'error8'] = abs(
            lampam_target[7]- result.lampam[7])
        table_result.loc[i, 'error9'] = abs(
            lampam_target[8]- result.lampam[8])
        table_result.loc[i, 'error10'] = abs(
            lampam_target[9]- result.lampam[9])
        table_result.loc[i, 'error11'] = abs(
            lampam_target[10]- result.lampam[10])
        table_result.loc[i, 'error12'] = abs(
            lampam_target[11]- result.lampam[11])

        # lampam_target
        table_result.loc[i, 'lampam_target[1]'] = lampam_target[0]
        table_result.loc[i, 'lampam_target[2]'] = lampam_target[1]
        table_result.loc[i, 'lampam_target[3]'] = lampam_target[2]
        table_result.loc[i, 'lampam_target[4]'] = lampam_target[3]
        table_result.loc[i, 'lampam_target[5]'] = lampam_target[4]
        table_result.loc[i, 'lampam_target[6]'] = lampam_target[5]
        table_result.loc[i, 'lampam_target[7]'] = lampam_target[6]
        table_result.loc[i, 'lampam_target[8]'] = lampam_target[7]
        table_result.loc[i, 'lampam_target[9]'] = lampam_target[8]
        table_result.loc[i, 'lampam_target[10]'] = lampam_target[9]
        table_result.loc[i, 'lampam_target[11]'] = lampam_target[10]
        table_result.loc[i, 'lampam_target[12]'] = lampam_target[11]

        # Retrieved stacking sequence at step 1
        ss_flatten = np.array(result.ss_tab[0], dtype=str)
        ss_flatten = ' '.join(ss_flatten)
        table_result.loc[i, 'ss retrieved at step 1'] = ss_flatten

        # Retrieved stacking sequence
        ss_flatten = np.array(result.ss, dtype=str)
        ss_flatten = ' '.join(ss_flatten)
        table_result.loc[i, 'ss retrieved'] = ss_flatten

        # Target stacking sequence
        ss_flatten = np.array(ss_target, dtype=str)
        #ss_flatten = ' '.join(ss_flatten)
        table_result.loc[i, 'ss target'] = ss_flatten

        # Ply counts
        table_result.loc[i, 'N0_target'] = N0_Target
        table_result.loc[i, 'N90_target'] = N90_Target
        table_result.loc[i, 'N45_target'] = N45_Target
        table_result.loc[i, 'N-45_target'] = N135_Target
        N0 = sum(result.ss == 0)
        N90 = sum(result.ss == 90)
        N45 = sum(result.ss == 45)
        N135 = sum(result.ss == -45)
        table_result.loc[i, 'N0 - N0_target'] = N0 - N0_Target
        table_result.loc[i, 'N90 - N90_target'] = N90 - N90_Target
        table_result.loc[i, 'N45 - N45_target'] = N45 - N45_Target
        table_result.loc[i, 'N-45 - N-45_target'] = N135 - N135_Target
        table_result.loc[i, 'penalty value for the 10% rule'] \
        = calc_penalty_10_ss(result.ss, constraints)

        for ind in range(n_outer_step):
            # numbers of stacks at the last level of the last group search
            table_result.loc[i, 'n_designs_last_level ' + str(ind + 1)] \
            = result.n_designs_last_level_tab[ind]
            # numbers of repaired stacks at the last group search
            table_result.loc[i, 'n_designs_repaired ' + str(ind + 1)] \
            = result.n_designs_repaired_tab[ind]
            # numbers of unique repaired stacks at the last group search
            table_result.loc[i, 'n_designs_repaired_unique ' + str(ind + 1)] \
            = result.n_designs_repaired_unique_tab[ind]

        # in-plane orthotropy
        ipo_now = ipo_param_1_12(result.lampam, mat_prop, constraints.sym)
        table_result.loc[i, 'In-plane orthotropy parameter 1'] = ipo_now[0]
        table_result.loc[i, 'In-plane orthotropy parameter 2'] = ipo_now[1]
        table_result.loc[i, 'In-plane orthotropy parameter 3'] = ipo_now[2]
        table_result.loc[i, 'In-plane orthotropy parameter 4'] = ipo_now[3]
        table_result.loc[i, 'In-plane orthotropy parameter 5'] = ipo_now[4]
        table_result.loc[i, 'In-plane orthotropy parameter 6'] = ipo_now[5]
        table_result.loc[i, 'In-plane orthotropy parameter 7'] = ipo_now[6]
        table_result.loc[i, 'In-plane orthotropy parameter 8'] = ipo_now[7]
        table_result.loc[i, 'In-plane orthotropy parameter 9'] = ipo_now[8]
        table_result.loc[i, 'In-plane orthotropy parameter 10'] = ipo_now[9]
        table_result.loc[i, 'In-plane orthotropy parameter 11'] = ipo_now[10]
        table_result.loc[i, 'In-plane orthotropy parameter 12'] = ipo_now[11]

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

#        print('A16', A16, A16_target, abs((A16 - A16_target)/A16_target))
#        print('A26', A26, A26_target, abs((A26 - A26_target)/A26_target))
#
#        print('D16', D16, D16_target, abs((D16 - D16_target)/D16_target))
#        print('D26', D26, D26_target, abs((D26 - D26_target)/D26_target))

        table_result.loc[i, 'diff A11 percentage'] \
        = abs((A11 - A11_target)/A11_target)
        table_result.loc[i, 'diff A22 percentage'] \
        = abs((A22 - A22_target)/A22_target)

        if abs(A12_target/A11_target) > 1e-8:
            table_result.loc[i, 'diff A12 percentage'] \
            = abs((A12 - A12_target)/A12_target)
        else:
            table_result.loc[i, 'diff A12 percentage'] = np.NaN
        if abs(A66_target/A11_target) > 1e-8:
            table_result.loc[i, 'diff A66 percentage'] \
            = abs((A66 - A66_target)/A66_target)
        else:
            table_result.loc[i, 'diff A66 percentage'] = np.NaN
        if abs(A16_target/A11_target) > 1e-8:
            table_result.loc[i, 'diff A16 percentage'] \
            = abs((A16 - A16_target)/A16_target)
        else:
            table_result.loc[i, 'diff A16 percentage'] = np.NaN
        if abs(A26_target/A11_target) > 1e-8:
            table_result.loc[i, 'diff A26 percentage'] \
            = abs((A26 - A26_target)/A26_target)
        else:
            table_result.loc[i, 'diff A26 percentage'] = np.NaN

        if B11_target:
            table_result.loc[i, 'diff B11 percentage'] \
            = abs((B11 - B11_target)/B11_target)
        else:
            table_result.loc[i, 'diff B11 percentage'] = np.NaN
        if B22_target:
            table_result.loc[i, 'diff B22 percentage'] \
            = abs((B22 - B22_target)/B22_target)
        else:
            table_result.loc[i, 'diff B22 percentage'] = np.NaN
        if B12_target:
            table_result.loc[i, 'diff B12 percentage'] \
            = abs((B12 - B12_target)/B12_target)
        else:
            table_result.loc[i, 'diff B12 percentage'] = np.NaN
        if B66_target:
            table_result.loc[i, 'diff B66 percentage'] \
            = abs((B66 - B66_target)/B66_target)
        else:
            table_result.loc[i, 'diff B66 percentage'] = np.NaN
        if B16_target:
            table_result.loc[i, 'diff B16 percentage'] \
            = abs((B16 - B16_target)/B16_target)
        else:
            table_result.loc[i, 'diff B16 percentage'] = np.NaN
        if B26_target:
            table_result.loc[i, 'diff B26 percentage'] \
            = abs((B26 - B26_target)/B26_target)
        else:
            table_result.loc[i, 'diff B26 percentage'] = np.NaN

        table_result.loc[i, 'diff D11 percentage'] \
        = abs((D11 - D11_target)/D11_target)
        table_result.loc[i, 'diff D22 percentage'] \
        = abs((D22 - D22_target)/D22_target)
        if abs(D12_target/D11_target) > 1e-8:
            table_result.loc[i, 'diff D12 percentage'] \
            = abs((D12 - D12_target)/D12_target)
        else:
            table_result.loc[i, 'diff D12 percentage'] = np.NaN
        if abs(D66_target/D11_target) > 1e-8:
            table_result.loc[i, 'diff D66 percentage'] \
            = abs((D66 - D66_target)/D66_target)
        else:
            table_result.loc[i, 'diff D66 percentage'] = np.NaN
        if abs(D16_target/D11_target) > 1e-8:
            table_result.loc[i, 'diff D16 percentage'] \
            = abs((D16 - D16_target)/D16_target)
        else:
            table_result.loc[i, 'diff D16 percentage'] = np.NaN
        if abs(D26_target/D11_target) > 1e-8:
            table_result.loc[i, 'diff D26 percentage'] \
            = abs((D26 - D26_target)/D26_target)
        else:
            table_result.loc[i, 'diff D26 percentage'] = np.NaN



#        table_result.loc[i, 'diff A11 percentage - approx'] \
#            = abs(4*(lampam_target[0] - result.lampam[0]) \
#                  + (lampam_target[1] - result.lampam[1])) \
#                  / abs(3 + 4*lampam_target[0] + lampam_target[1])
#        table_result.loc[i, 'diff A22 percentage - approx'] \
#            = abs(4*(lampam_target[0] - result.lampam[0]) \
#                  - (lampam_target[1] - result.lampam[1])) \
#                  / abs(3 - 4*lampam_target[0] + lampam_target[1])
#        print(lampam_target[0], lampam_target[1])
#        if abs(A12_target/A11_target) > 1e-8:
#            table_result.loc[i, 'diff A12 percentage - approx'] \
#            = abs((lampam_target[1] - result.lampam[1])) \
#                  / abs(1 - lampam_target[1])
#        else:
#            table_result.loc[i, 'diff A12 percentage - approx'] = np.NaN
#        if abs(A66_target/A11_target) > 1e-8:
#            table_result.loc[i, 'diff A66 percentage - approx'] \
#            = abs((lampam_target[1] - result.lampam[1])) \
#                  / abs(4 - lampam_target[1])
#        else:
#            table_result.loc[i, 'diff A66 percentage - approx'] = np.NaN
#        if abs(A16_target/A11_target) > 1e-8:
#            table_result.loc[i, 'diff A16 percentage - approx'] \
#            = abs(2*(lampam_target[2] - result.lampam[2]) \
#                  + (lampam_target[3] - result.lampam[3])) \
#                  / abs(2*lampam_target[2] + lampam_target[3])
#        else:
#            table_result.loc[i, 'diff A16 percentage - approx'] = np.NaN
#        if abs(A26_target/A11_target) > 1e-8:
#            table_result.loc[i, 'diff A26 percentage - approx'] \
#            = abs(2*(lampam_target[2] - result.lampam[2]) \
#                  - (lampam_target[3] - result.lampam[3])) \
#                  / abs(2*lampam_target[2] - lampam_target[3])
#        else:
#            table_result.loc[i, 'diff A26 percentage - approx'] = np.NaN
#
#
#
#        table_result.loc[i, 'diff D11 percentage - approx'] \
#            = abs(4*(lampam_target[8] - result.lampam[8]) \
#                  + (lampam_target[9] - result.lampam[9])) \
#                  / abs(3 + 4*lampam_target[8] + lampam_target[9])
#        table_result.loc[i, 'diff D22 percentage - approx'] \
#            = abs(4*(lampam_target[8] - result.lampam[8]) \
#                  - (lampam_target[9] - result.lampam[9])) \
#                  / abs(3 - 4*lampam_target[8] + lampam_target[9])
#        if abs(D12_target/D11_target) > 1e-8:
#            table_result.loc[i, 'diff D12 percentage - approx'] \
#            = abs((lampam_target[9] - result.lampam[9])) \
#                  / abs(1 - lampam_target[9])
#        else:
#            table_result.loc[i, 'diff D12 percentage - approx'] = np.NaN
#        if abs(D66_target/D11_target) > 1e-8:
#            table_result.loc[i, 'diff D66 percentage - approx'] \
#            = abs((lampam_target[9] - result.lampam[9])) \
#                  / abs(4 - lampam_target[9])
#        else:
#            table_result.loc[i, 'diff D66 percentage - approx'] = np.NaN
#        if abs(D16_target/D11_target) > 1e-8:
#            table_result.loc[i, 'diff D16 percentage - approx'] \
#            = abs(2*(lampam_target[10] - result.lampam[10]) \
#                  + (lampam_target[11] - result.lampam[11])) \
#                  / abs(2*lampam_target[10] + lampam_target[11])
#        else:
#            table_result.loc[i, 'diff D16 percentage - approx'] = np.NaN
#        if abs(D26_target/D11_target) > 1e-8:
#            table_result.loc[i, 'diff D26 percentage - approx'] \
#            = abs(2*(lampam_target[10] - result.lampam[10]) \
#                  - (lampam_target[11] - result.lampam[11])) \
#                  / abs(2*lampam_target[10] - lampam_target[11])
#        else:
#            table_result.loc[i, 'diff D26 percentage - approx'] = np.NaN


### Write results in a excell sheet
writer = pd.ExcelWriter(result_filename)
table_result.to_excel(writer, 'results')
writer.save()
save_constraints_LAYLA(result_filename, constraints)
save_parameters_LAYLA_V02(result_filename, parameters)
save_materials(result_filename, mat_prop)
autofit_column_widths(result_filename)
