# -*- coding: utf-8 -*-
"""
This script creates the iput files used for testing BELLA
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
from src.BELLA.format_pdl import convert_sst_to_ss
from src.BELLA.save_set_up import save_constraints_BELLA
from src.BELLA.save_set_up import save_multipanel
from src.BELLA.save_set_up import save_objective_function_BELLA
from src.BELLA.save_set_up import save_materials
from src.guidelines.one_stack import check_lay_up_rules
from src.guidelines.one_stack import check_ply_drop_rules
from src.divers.excel import delete_file, autofit_column_widths

sheet = 'SST-40-60'
sheet = 'SST-80-120'
sheet = 'SST-120-180'

filename_input = '/BELLA/input-files/SST.xlsx'
filename_res = 'input_file_' + sheet + '.xlsx'

# check for authorisation before overwriting
delete_file(filename_res)

# number of panels
if sheet == 'SST-40-60': n_panels = 6
elif sheet == 'SST-80-120': n_panels = 11
elif sheet == 'SST-120-180': n_panels = 16

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
if sheet == 'SST-40-60':
    # panel IDs
    ID = np.arange(1, 7)
    # number of plies in each panel
    n_plies_per_panel = np.arange(40, 61, 4)
elif sheet == 'SST-80-120':
    # panel IDs
    ID = np.arange(1, 12)
    # number of plies in each panel
    n_plies_per_panel = np.arange(80, 121, 4)
elif sheet == 'SST-120-180':
    # panel IDs
    ID = np.arange(1, 17)
    # number of plies in each panel
    n_plies_per_panel = np.arange(120, 181, 4)

# panels adjacency
neighbour_panels = dict()
neighbour_panels[1] = [2]
neighbour_panels[n_panels] = [n_panels - 1]
for i in range(2, n_panels):
    neighbour_panels[i] = [i - 1, i + 1]

sst = pd.read_excel(filename_input, sheet_name=sheet).fillna(-1)
sst = np.array(sst, int).T
sst = np.hstack((sst, np.flip(sst, axis=1)))
sst = sst[::2]

ss = convert_sst_to_ss(sst)

# panel amination parameters targets
lampam_targets = [calc_lampam(stack) for stack in ss]

panels = []
for ind_panel in range(n_panels):
    panels.append(Panel(
        ID=ID[ind_panel],
        lampam_target=lampam_targets[ind_panel],
        lampam_weightings=lampam_weightings,
        n_plies=len(ss[ind_panel]),
        length_x=0,
        length_y=0,
        N_x=0,
        N_y=0,
        area=0,
        weighting=panel_weightings[ind_panel],
        neighbour_panels=neighbour_panels[ID[ind_panel]],
        constraints=constraints))

multipanel = MultiPanel(panels)
multipanel.filter_target_lampams(constraints, obj_func_param)
multipanel.filter_lampam_weightings(constraints, obj_func_param)

### Checks for feasibility of the multi-panel composite layout ----------------

check_ply_drop_rules(sst, multipanel, constraints, reduced=False)

for stack in ss:
    check_lay_up_rules(stack, constraints)

### Save data -----------------------------------------------------------------

save_multipanel(filename_res, multipanel, obj_func_param, calc_penalties=True,
                constraints=constraints, sst=sst)
save_constraints_BELLA(filename_res, constraints)
save_objective_function_BELLA(filename_res, obj_func_param)
save_materials(filename_res, materials)
autofit_column_widths(filename_res)