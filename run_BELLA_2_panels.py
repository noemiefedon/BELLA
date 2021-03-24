# -*- coding: utf-8 -*-
"""
Script to retrieve a blended multi-panel layout based on:
    - a panel thickness distribution
    - set of lamination parameter targets for each panel
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import time
import numpy as np
import random

random.seed(10)

sys.path.append(r'C:\BELLA')
from src.BELLA.panels import Panel
from src.BELLA.multipanels import MultiPanel
from src.BELLA.parameters import Parameters
from src.BELLA.constraints import Constraints
from src.BELLA.obj_function import ObjFunction
from src.BELLA.materials import Material
from src.BELLA.optimiser import BELLA_optimiser
from src.BELLA.format_pdl import convert_sst_to_ss
from src.CLA.lampam_functions import calc_lampam_2
from src.divers.excel import delete_file
from src.divers.pretty_print import print_list_ss

filename = 'BELLA-2-panels-test.xlsx'

# check for authorisation before overwriting
delete_file(filename)

constraints_set = 'C0'
constraints_set = 'C1'

### Targets and panel geometries ----------------------------------------------

# panel IDs
ID = [1, 2]

# number of panels
n_panels = len(ID)

# panels adjacency
neighbour_panels = {1:[2], 2:[1]}

# boundary weights
boundary_weights = {(1, 2) : 1}

# panel length in the x-direction
length_x = [1] *n_panels

# panel length in the y-direction
length_y = [1] *n_panels

# panel loading per unit width in the x-direction
N_x = [1] *n_panels

# panel loading per unit width in the y-direction
N_y = [1] *n_panels

# panel target stacking sequences
sst_target = np.array([
    45, 0, -45, 0, 0, 45, 0, 0, 45, 90, 90, -45, -45, 0, ], dtype=int)
sst_target = np.array([
    0, 0, 0, 0, 0, 0, 0, 0, ], dtype=int)
sst_target = np.hstack((sst_target, np.flip(sst_target)))
sst_target = np.vstack((sst_target, np.copy(sst_target), np.copy(sst_target)))
sst_target[1:, 3] = -1
sst_target[1:, 6] = -1
# sst_target[1:, 9] = -1
sst_target[1:, -4] = -1
sst_target[1:, -7] = -1
# sst_target[1:, -10] = -1
# sst_target[2:, 4] = -1
# sst_target[2:, 10] = -1
# sst_target[2:, -5] = -1
# sst_target[2:, -11] = -1

ss_target = convert_sst_to_ss(sst_target)

# panel number of plies
n_plies = [stack.size for stack in ss_target]

# panel amination parameters targets
lampam_targets = [
    calc_lampam_2(ss_target[ind_panel]) for ind_panel in range(n_panels)]

### Design guidelines ---------------------------------------------------------

# constraints_set == 'C0' ->
#   - ply-drop spacing rule enforced with a minimum of
#   constraints.min_drop plies between ply drops at panel boundaries
#   - covering rule enforced by preventing the drop of the
#    constraints.n_covering outermost plies on each laminate surface
#   - symmetry rule enforced, no other lay-up rules
#
# constraints_set == 'C1' ->
#   - ply-drop spacing rule enforced with a minimum of
#   constraints.min_drop plies between ply drops at panel boundaries
#   - covering enforrced by preventing the drop of the
#    constraints.n_covering outermost plies on each laminate surface
#   - symmetry rule enforced
#   - 10% rule enforced
#      if rule_10_Abdalla == True rule applied by restricting LPs instead of
#          ply percentages and percent_Abdalla is the percentage limit of the
#          rule
#      otherwise:
#         if combined_45_135 == True the restrictions are:
#           - a maximum percentage of constraints.percent_0 0 deg plies
#           - a maximum percentage of constraints.percent_90 90 deg plies
#           - a maximum percentage of constraints.percent_45_135 +-45 deg plies
#         if combined_45_135 == False the restrictions are:
#           - a maximum percentage of constraints.percent_0 0 deg plies
#           - a maximum percentage of constraints.percent_90 90 deg plies
#           - a maximum percentage of constraints.percent_45 45 deg plies
#           - a maximum percentage of constraints.percent_135 -45 deg plies
#   - disorientation rule enforced with variation of fibre angle between
#    adacent plies limited to a maximum value of constraints.delta_angle
#    degrees
#   - contiguity rule enforced with no more than constraints.n_contig
#    adajacent plies with same fibre angle
#   - damage tolerance rule enforced
#     if constraints.dam_tol_rule == 1 the restrictions are:
#       - one outer ply at + or -45 deg at the laminate surfaces
#        (2 plies intotal)
#     if constraints.dam_tol_rule == 2 the restrictions are:
#       - [+45, -45] or [-45, +45] at the laminate surfaces
#        (4 plies in total)
#     if constraints.dam_tol_rule == 3 the restrictions are:
#       - [+45,-45] [-45,+45] [+45,+45] or [-45,-45] at the laminate
#         surfaces (4 plies in total)
#   - out-of-plane orthotropy rule enforced to have small absolutes values
#    of LP_11 and LP_12 such that the values of D16 and D26 are small too

## lay-up rules

# set of admissible fibre orientations
set_of_angles = np.array([-45, 0, 45, 90], dtype=int)
#set_of_angles = np.array([-45, 0, 45, 90, +30, -30, +60, -60], dtype=int)

sym = True # symmetry rule
oopo = False # out-of-plane orthotropy requirements

if constraints_set == 'C0':
    bal = False # balance rule
    rule_10_percent = False # 10% rule
    diso = False # disorientation rule
    contig = False # contiguity rule
    dam_tol = False # damage-tolerance rule
else:
    bal = True
    rule_10_percent = True
    diso = True
    contig = True
    dam_tol = True

rule_10_Abdalla = False # 10% rule restricting LPs instead of ply percentages
percent_Abdalla = 0 # percentage limit for the 10% rule applied on LPs
combine_45_135 = True # True if restriction on +-45 plies combined for 10% rule
percent_0 = 10 # percentage used in the 10% rule for 0 deg plies
percent_45 = 0 # percentage used in the 10% rule for +45 deg plies
percent_90 = 10 # percentage used in the 10% rule for 90 deg plies
percent_135 = 0 # percentage used in the 10% rule for -45 deg plies
percent_45_135 =10 # percentage used in the 10% rule for +-45 deg plies
delta_angle = 45 # maximum angle difference for adjacent plies
n_contig = 5 # maximum number of adjacent plies with the same fibre orientation
dam_tol_rule = 1 # type of damage tolerance rule

## ply-drop rules

covering = True # covering rule
n_covering = 1 # number of plies ruled by covering rule at laminate surfaces
pdl_spacing = True # ply drop spacing rule
min_drop = 2 # Minimum number of continuous plies between ply drops

constraints = Constraints(
    sym=sym,
    bal=bal,
    oopo=oopo,
    dam_tol=dam_tol,
    dam_tol_rule=dam_tol_rule,
    covering=covering,
    n_covering=n_covering,
    rule_10_percent=rule_10_percent,
    rule_10_Abdalla=rule_10_Abdalla,
    percent_Abdalla=percent_Abdalla,
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
    set_of_angles=set_of_angles,
    min_drop=min_drop,
    pdl_spacing=pdl_spacing)

### Material properties -------------------------------------------------------

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

### Objective function parameters ---------------------------------------------

# Coefficient for the 10% rule penalty
coeff_10 = 1
# Coefficient for the contiguity constraint penalty
coeff_contig = 1
# Coefficient for the disorientation constraint penalty
coeff_diso = 1
# Coefficient for the out-of-plane orthotropy penalty
coeff_oopo = 1
# Coefficient for the ply drop spacing guideline penalty
coeff_spacing = 1

# Lamination-parameter weightings in panel objective functions
# (In practice these weightings can be different for each panel)
optimisation_type = 'AD'
if optimisation_type == 'A':
    if all(elem in {0, +45, -45, 90}  for elem in constraints.set_of_angles):
        lampam_weightings = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    else:
        lampam_weightings = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
elif optimisation_type == 'D':
    if all(elem in {0, +45, -45, 90}  for elem in constraints.set_of_angles):
        lampam_weightings = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0])
    else:
        lampam_weightings = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1])
elif optimisation_type == 'AD':
    if all(elem in {0, +45, -45, 90}  for elem in constraints.set_of_angles):
        lampam_weightings = np.array([1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0])
    else:
        lampam_weightings = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])

# Weightings of the panels in the multi-panel objecive function
panel_weightings = np.ones((n_panels,), float)

obj_func_param = ObjFunction(
    constraints=constraints,
    coeff_contig=coeff_contig,
    coeff_diso=coeff_diso,
    coeff_10=coeff_10,
    coeff_oopo=coeff_oopo,
    coeff_spacing=coeff_spacing)

### Optimiser  ----------------------------------------------------------------

## For step 2 of BELLA

# Number of initial ply drops to be tested
n_ini_ply_drops = 1
# Minimum ply count for ply groups during ply drop layout generation
group_size_min = 4
# Desired ply count for ply groups during ply drop layout generation
group_size_max = 12
# Time limit to create a group ply-drop layout
time_limit_group_pdl = 1
# Time limit to create a ply-drop layout
time_limit_all_pdls = 10

## For step 3 of BELLA

# Branching limit for global pruning
global_node_limit = 3
# Branching limit for local pruning
local_node_limit = 3
# Branching limit for global pruning at the last level 
global_node_limit_final = 1
# Branching limit for local pruning at the last level 
local_node_limit_final = 1

## For step 4.1 of BELLA

# to save repair success rates
save_success_rate = True
# Thickness of the reference panels
n_plies_ref_panel = 1
# repair to improve the convergence towards the in-plane lamination parameter
# targets
repair_membrane_switch = False
# repair to improve the convergence towards the out-of-plane lamination
# parameter targets
repair_flexural_switch = False
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

## For step 4.2 of BELLA

# Branching limit for global pruning during ply drop layout optimisation
global_node_limit2 = 1
# Branching limit for local pruning during ply drop layout optimisation
local_node_limit2 = 1

## For step 4.3 of BELLA

# Branching limit for global pruning during ply drop layout optimisation
global_node_limit3 = 1
# Branching limit for local pruning during ply drop layout optimisation
local_node_limit3 = 1

parameters = Parameters(
    constraints=constraints,
    group_size_min=group_size_min,
    group_size_max=group_size_max,
    n_ini_ply_drops=n_ini_ply_drops,
    global_node_limit=global_node_limit,
    global_node_limit_final=global_node_limit_final,
    local_node_limit=local_node_limit,
    local_node_limit_final=local_node_limit_final,
    global_node_limit2=global_node_limit2,
    local_node_limit2=local_node_limit2,
    global_node_limit3=global_node_limit3,
    local_node_limit3=local_node_limit3,
    save_success_rate=save_success_rate,
    p_A=p_A,
    n_D1=n_D1,
    n_D2=n_D2,
    n_D3=n_D3,
    repair_membrane_switch=repair_membrane_switch,
    repair_flexural_switch=repair_flexural_switch,
    n_plies_ref_panel=n_plies_ref_panel,
    time_limit_group_pdl=time_limit_group_pdl,
    time_limit_all_pdls=time_limit_all_pdls)

### Multi-panel composite laminate layout -------------------------------------

panels = []
for ind_panel in range(n_panels):
    panels.append(Panel(
        ID=ID[ind_panel],
        lampam_target=lampam_targets[ind_panel],
        lampam_weightings=lampam_weightings,
        n_plies=n_plies[ind_panel],
        length_x=length_x[ind_panel],
        length_y=length_y[ind_panel],
        N_x=N_x[ind_panel],
        N_y=N_y[ind_panel],
        weighting=panel_weightings[ind_panel],
        neighbour_panels=neighbour_panels[ID[ind_panel]],
        constraints=constraints))


multipanel = MultiPanel(panels, boundary_weights)

### Optimiser Run -------------------------------------------------------------

t = time.time()
result = BELLA_optimiser(multipanel, parameters, obj_func_param, constraints,
                         filename, mat)
elapsed = time.time() - t

### Display results -----------------------------------------------------------

# solution stacking sequences (list of arrays)
ss = result.ss
# solution stacking sequence table (array)
sst = result.sst
# solution lamination parameters (array)
lampam = result.lampam
# objective function value of the solution with constraint
# penalties
obj_constraints = result.obj_constraints
# objective function value of the solution with no constraint
# penalties
obj_no_constraints = result.obj_no_constraints
# objective function values with no constraint
#penalties
obj_no_constraints_tab = result.obj_no_constraints_tab
# objective function values with constraint
# penalties
obj_constraints_tab = result.obj_constraints_tab
# number of calls of the objective functions
n_obj_func_calls_tab = result.n_obj_func_calls_tab
# penalty for the disorientation constraint
penalty_diso_tab = result.penalty_diso_tab
# penalty for the contiguity constraint
penalty_contig_tab = result.penalty_contig_tab
# penalty for the 10% rule
penalty_10_tab = result.penalty_10_tab
# penalty for in-plane orthotropy
penalty_bal_ipo_tab = result.penalty_bal_ipo_tab
# penalty for out-of-plane orthotropy
penalty_oopo_tab = result.penalty_oopo_tab
# index for the outer loop with the best design
ind_mini = result.ind_mini
# number of plies in each fibre direction for each panel
n_plies_per_angle = result.n_plies_per_angle

print(r'\\\\\\\ Final objective : ', obj_constraints)
print(r'\\\\\\\ Elapsed time : ', elapsed, 's')
print(r'\\\\\\\ objectives: ', obj_constraints_tab)
print(r'\\\\\\\ Number of function evaluations')
print(n_obj_func_calls_tab)
if constraints.rule_10_percent:
    print(r'\\\\\\\ Penalties for the 10% rule')
    print(penalty_10_tab)
if constraints.diso:
    print(r'\\\\\\\ Penalties for disorientation')
    print(penalty_diso_tab)
if constraints.contig:
    print(r'\\\\\\\ Penalties for contiguity')
    print(penalty_contig_tab)
if constraints.ipo:
    print(r'\\\\\\\ Penalties for in-plane-orthotropy')
    print(penalty_bal_ipo_tab)
if constraints.oopo:
    print(r'\\\\\\\ Penalties for out-of-plane-orthotropy')
    print(penalty_oopo_tab)

print('\nRetrieved stacking sequences')
#    for ii, panel in enumerate(multipanel.panels):
#        print('panel', ii + 1)
#        print_ss(ss[ii])

#print('lampam_Retrieved VS lampam_target & difference')
#for ii, panel in enumerate(multipanel.panels):
#    print('panel', ii + 1)
#    print_lampam(lampam[ii], panel.lampam_target, diff=True)

#for index, panel in enumerate(multipanel.panels):
#    print('panel', index + 1)
#    N15plus = sum(ss[index] == 15)
#    N15minus = sum(ss[index] == -15)
#    N30plus = sum(ss[index] == 30)
#    N30minus = sum(ss[index] == -30)
#    N45plus = sum(ss[index] == 45)
#    N45minus = sum(ss[index] == -45)
#    N60plus = sum(ss[index] == 60)
#    N60minus = sum(ss[index] == -60)
#    N75plus = sum(ss[index] == 15)
#    N75minus = sum(ss[index] == -15)
    # balance check
#    if N15plus != N15minus \
#    or N30plus != N30minus \
#    or N45plus != N45minus \
#    or N60plus != N60minus \
#    or N75plus != N75minus:
#    if N45plus != N45minus:
#        print('balance constraint not respected:')
#        print(f'N45 - N-45: {N45plus - N45minus}')

print_list_ss(ss)