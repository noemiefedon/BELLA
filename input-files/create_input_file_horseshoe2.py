# -*- coding: utf-8 -*-
"""
This script saves the input file for the feasible horseshoe problem.
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'


import sys
import pandas as pd
import numpy as np
import numpy.matlib

sys.path.append(r'C:\BELLA')
from src.CLA.lampam_functions import calc_lampam
from src.BELLA.panels import Panel
from src.BELLA.multipanels import MultiPanel
from src.BELLA.constraints import Constraints
from src.BELLA.obj_function import ObjFunction
from src.BELLA.materials import Material
from src.BELLA.save_set_up import save_constraints_BELLA
from src.BELLA.save_set_up import save_multipanel
from src.BELLA.save_set_up import save_objective_function_BELLA
from src.BELLA.save_set_up import save_materials
from src.guidelines.one_stack import check_lay_up_rules
from src.guidelines.one_stack import check_ply_drop_rules
from src.divers.excel import delete_file, autofit_column_widths

sheet = 'horseshoe2'
filename_input = '/BELLA/input-files/SST.xlsx'
filename_res = 'input_file_' + sheet + '.xlsx'

# check for authorisation before overwriting
delete_file(filename_res)

# number of panels
n_panels = 18

### Design guidelines ---------------------------------------------------------

constraints_set = 'C0'
constraints_set = 'C1'

## lay-up rules

# set of admissible fibre orientations
set_of_angles = np.array([-45, 0, 45, 90], dtype=int)
set_of_angles = np.array([
    -45, 0, 45, 90, +30, -30, +60, -60, 15, -15, 75, -75], dtype=int)

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

rule_10_Abdalla = True # 10% rule restricting LPs instead of ply percentages
percent_Abdalla = 10 # percentage limit for the 10% rule applied on LPs
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

### Objective function parameters ---------------------------------------------

# Coefficient for the 10% rule penalty
coeff_10 = 1
# Coefficient for the contiguity constraint penalty
coeff_contig = 1
# Coefficient for the disorientation constraint penalty
coeff_diso = 10
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

## Multi-panel objective function

# Weightings of the panels in the multi-panel objecive function
panel_weightings = np.ones((n_panels,), float)

obj_func_param = ObjFunction(
    constraints=constraints,
    coeff_contig=coeff_contig,
    coeff_diso=coeff_diso,
    coeff_10=coeff_10,
    coeff_oopo=coeff_oopo,
    coeff_spacing=coeff_spacing)

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

materials = Material(E11=E11, E22=E22, G12=G12, nu12=nu12,
               density_area=density_area, ply_t=ply_t)

### Multi-panel composite laminate layout -------------------------------------

# panel IDs
ID = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]

# panel number of plies
n_plies_per_panel = [32, 28, 20, 18, 16, 22, 18, 24, 38,
                     34, 30, 28, 22, 18, 24, 30, 18, 22]
# panels adjacency
neighbour_panels = {
    1 : [2, 9],
    2 : [1, 3, 6, 10],
    3 : [2, 4, 6],
    4 : [3, 5, 7],
    5 : [4, 8],
    6 : [2, 3, 7],
    7 : [4, 6, 8],
    8 : [5, 7],
    9 : [1, 10, 11],
    10 : [2, 9, 12],
    11 : [9, 12],
    12 : [10, 11, 13, 16],
    13 : [12, 14, 16],
    14 : [13, 15, 17],
    15 : [14, 18],
    16 : [12, 13, 17],
    17 : [14, 16, 18],
    18 : [15, 17]}

# boundary weights
boundary_weights = {(1, 2) : 0.610,
                    (1, 9) : 0.457,
                    (2, 3) : 0.305,
                    (2, 6) : 0.305,
                    (2, 10) : 0.457,
                    (3, 4) : 0.305,
                    (3, 6) : 0.508,
                    (4, 5) : 0.305,
                    (4, 7) : 0.508,
                    (5, 8) : 0.508,
                    (6, 7) : 0.305,
                    (7, 8) : 0.305,
                    (9, 10) : 0.610,
                    (9, 11) : 0.457,
                    (10, 12) : 0.457,
                    (11, 12) : 0.610,
                    (12, 13) : 0.305,
                    (12, 16) : 0.305,
                    (13, 14) : 0.305,
                    (13, 16) : 0.508,
                    (14, 15) : 0.305,
                    (14, 17) : 0.508,
                    (15, 18) : 0.508,
                    (16, 17) : 0.305,
                    (17, 18) : 0.305}

# panel length in the x-direction (m)
length_x = (25.40/1000)*np.array([18, 18, 20, 20, 20, 20, 20, 20,
    18, 18, 18, 18, 20, 20, 20, 20, 20, 20])

# panel length in the y-direction (m)
length_y = (25.40/1000)*np.array([24, 24, 12, 12, 12, 12, 12, 12,
    24, 24, 24, 24, 12, 12, 12, 12, 12, 12])

# 1 lbf/in = 0.175127 N/mm
# panel loading per unit width in the x-direction in N/m
N_x = 175.127*np.array([700, 375, 270, 250, 210, 305, 290, 600,
                        1100, 900, 375, 400, 330, 190, 300, 815, 320, 300])

# panel loading per unit width in the y-direction in N/m
N_y = 175.127*np.array([400, 360, 325, 200, 100, 360, 195, 480,
                        600, 400, 525, 320, 330, 205, 610, 1000, 180, 410])


sst = pd.read_excel(filename_input, sheet_name=sheet).fillna(-1)
sst = np.array(sst, int).T
sst = np.hstack((sst, np.flip(sst, axis=1)))

n_plies_2_ss = dict()
n_plies_2_lampam = dict()
n_plies_2_sst = dict()

for stack_sst in sst:
    stack = np.copy(stack_sst)
    for ind_ply in range(38)[::-1]:
        if stack[ind_ply] == -1:
            stack = np.delete(stack, ind_ply)
    n_plies_2_sst[len(stack)] = stack_sst
    n_plies_2_ss[len(stack)] = stack
    n_plies_2_lampam[len(stack)] = calc_lampam(stack)

sst = np.array([n_plies_2_sst[n] for n in n_plies_per_panel])

panels = []
for ind_panel in range(n_panels):
    panels.append(Panel(
        ID=ID[ind_panel],
        lampam_target=n_plies_2_lampam[n_plies_per_panel[ind_panel]],
        lampam_weightings=lampam_weightings,
        n_plies=n_plies_per_panel[ind_panel],
        length_x=length_x[ind_panel],
        length_y=length_y[ind_panel],
        N_x=N_x[ind_panel],
        N_y=N_y[ind_panel],
        weighting=panel_weightings[ind_panel],
        neighbour_panels=neighbour_panels[ID[ind_panel]],
        constraints=constraints))

multipanel = MultiPanel(panels)
multipanel.filter_target_lampams(constraints, obj_func_param)
multipanel.filter_lampam_weightings(constraints, obj_func_param)

### Checks for feasibility of the multi-panel composite layout ----------------

#check_ply_drop_rules(sst, multipanel, constraints, reduced=False)

for stack in n_plies_2_ss.values():
    check_lay_up_rules(stack, constraints,
                       no_ipo_check=True, no_bal_check=True)

### Save data -----------------------------------------------------------------

save_multipanel(filename_res, multipanel, obj_func_param, calc_penalties=True,
                constraints=constraints, sst=sst)
save_constraints_BELLA(filename_res, constraints)
save_objective_function_BELLA(filename_res, obj_func_param)
save_materials(filename_res, materials)
autofit_column_widths(filename_res)

