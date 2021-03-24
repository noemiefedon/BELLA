# -*- coding: utf-8 -*-
"""
This script retrieves LAminate LAY-ups from one set of lamination parameters.
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import time

import numpy as np
import numpy.matlib
import math as ma

sys.path.append(r'C:\BELLA_and_LAYLA')
from src.LAYLA_V02.targets import Targets
from src.LAYLA_V02.parameters import Parameters
from src.LAYLA_V02.constraints import Constraints
from src.LAYLA_V02.optimiser import LAYLA_optimiser
from src.CLA.lampam_functions import calc_lampam_2
from src.BELLA.materials import Material
from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

#==============================================================================
# Targets
#==============================================================================
# Total number of plies
n_plies = 61
# Stacking sequence target
ss_target = np.array([45, 90, -45, 0, 0], dtype='int16')
print_ss(ss_target, 200)

# Calculation of target lamination parameters
lampam = calc_lampam_2(ss_target)
print_lampam(lampam)

targets = Targets(n_plies=n_plies, lampam=lampam, stack=ss_target)

#==============================================================================
# Type of optimisations
#==============================================================================
optimisation_type = 'A' # only in-plane lamination parameters optimised
optimisation_type = 'D' # only out-of-plane lamination parameters optimised
optimisation_type = 'AD' # in- and out-of-plane lamination parameters optimised

#==============================================================================
# Design guidelines
#==============================================================================
### Set of design and manufacturing constraints:
constraints_set = 'C0'
constraints_set = 'C1'
#     C0:   - No design and manufacturing constraints other than symmetry
#     C1:   - in-plane orthotropy enforced with penalties and repair
#           - 10% rule enforced with repair
#                   - 10% 0deg plies
#                   - 10% 90 deg plies
#                   - 10% 45deg plies
#                   - 10% -45 deg plies
#           - disorientation rule with Delta(theta) = 45 deg
#           - contiguity rule with n_contig = 5

# set of admissible fibre orientations
set_of_angles = np.array([-45, 0, 45, 90], dtype=int)
set_of_angles = np.array([-45, 0, 45, 90, +30, -30, +60, -60], dtype=int)

# symmetry
sym = False

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
dam_tol = True
# rule 1: one outer ply at + or -45 deg at laminate surfaces
# rule 2: [+45, -45] or [-45, +45] plies at laminate surfaces
# rule 3: [+45, -45], [+45, +45], [-45, -45] or [-45, +45] plies at laminate
dam_tol_rule = 1
dam_tol_rule = 2
#dam_tol_rule = 3

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

n_contig = 4
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
n_outer_step = 1

# branching limit for global pruning during ply orientation optimisation
global_node_limit = 8
# branching limit for local pruning during ply orientation optimisation
local_node_limit = 8
# branching limit for global pruning at the penultimate level during ply
# orientation optimisation
global_node_limit_p = 8
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
penalty_10_lampam_switch = True
# penalty for in-plane orthotropy, based on lamination parameters
penalty_ipo_switch = True
# penalty for balance, based on ply counts
penalty_bal_switch = False
# balanced laminate scheme
balanced_scheme = False

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
group_size_max = np.array([1000, 12, 12, 12, 12])

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
#     Optimiser Run
#==============================================================================
print('Algorithm running')

#print(targets)
#print(mat_prop)
#print(constraints)
#print(parameters)

t = time.time()
result = LAYLA_optimiser(parameters, constraints, targets, mat_prop)
elapsed1 = time.time() - t

#==============================================================================
#     Results Display
#==============================================================================
print()

print('\\\\\\\\\\\ objective with modified lamination parameters: ',
      result.objective)

print('\\\\\\\\\\\ Elapsed time : ', elapsed1, 's')

print('\nRetrieved stacking sequence')
print_ss(result.ss)

print('Difference of lamination parameters')
print_lampam(result.lampam - targets.lampam)

print('\nRetrieved lamination parameters')
print_lampam(result.lampam)

print(f'\nNumber of outer loops performed: {result.number_of_outer_steps_performed}')

print('result.n_designs_repaired_unique_tab',
      result.n_designs_repaired_unique_tab)