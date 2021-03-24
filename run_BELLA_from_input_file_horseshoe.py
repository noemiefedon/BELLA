# -*- coding: utf-8 -*-
"""
Script to retrieve a blended multi-panel layout based on:
    - a panel thickness distribution
    - set of lamination parameter targets for each panel
The thicknesses and lamination parameters match the targets of the horseshoe
problem and are found in input_file_horseshoe_... .xlsx.
y"""

__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import pandas as pd
import numpy as np
import numpy.matlib
import random
random.seed(0)

sys.path.append(r'C:\BELLA')
from src.BELLA.panels import Panel
from src.BELLA.multipanels import MultiPanel
from src.BELLA.parameters import Parameters
from src.BELLA.constraints import Constraints
from src.BELLA.obj_function import ObjFunction
from src.BELLA.materials import Material
from src.BELLA.optimiser import BELLA_optimiser
from src.divers.excel import delete_file

filename = 'input_file_horseshoe.xlsx'
filename = 'input_file_horseshoe2.xlsx'
filename_input = '/BELLA/input-files/' + filename
filename_result = ('results_BELLA_thin_' + filename).replace('input_file_', '')

# check for authorisation before overwriting
delete_file(filename_result)

### Design guidelines ---------------------------------------------------------

data_constraints = pd.read_excel(filename_input, sheet_name='Constraints',
                                 header=None, index_col=0).T
sym = data_constraints["symmetry"].iloc[0]
bal = data_constraints["balance"].iloc[0]
oopo = data_constraints["out-of-plane orthotropy"].iloc[0]
dam_tol = data_constraints["damage tolerance"].iloc[0]
dam_tol_rule = int(data_constraints["dam_tol_rule"].iloc[0])
covering = data_constraints["covering"].iloc[0]
n_covering = int(data_constraints["n_covering"].iloc[0])
rule_10_percent = data_constraints["10% rule"].iloc[0]
rule_10_Abdalla = data_constraints["10% rule applied on LPs"].iloc[0]
percent_Abdalla = float(data_constraints[
    "percentage limit when rule applied on LPs"].iloc[0])
percent_0 = float(data_constraints["percent_0"].iloc[0])
percent_45 = float(data_constraints["percent_45"].iloc[0])
percent_90 = float(data_constraints["percent_90"].iloc[0])
percent_135 = float(data_constraints["percent_-45"].iloc[0])
percent_45_135 = float(data_constraints["percent_+-45"].iloc[0])
diso = data_constraints["diso"].iloc[0]
delta_angle = float(data_constraints["delta_angle"].iloc[0])
contig = data_constraints["contig"].iloc[0]
n_contig = int(data_constraints["n_contig"].iloc[0])
set_of_angles = np.array(
    data_constraints["fibre orientations"].iloc[0].split(" "), int)
pdl_spacing = data_constraints["ply drop spacing rule"].iloc[0]
min_drop = int(data_constraints[
    "minimum number of continuous plies between ply drops"].iloc[0])
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
    diso=diso,
    contig=contig,
    n_contig=n_contig,
    delta_angle=delta_angle,
    set_of_angles=set_of_angles,
    min_drop=min_drop,
    pdl_spacing=pdl_spacing)

### Optimiser parameters ------------------------------------------------------

## For step 2 of BELLA

# Number of initial ply drops to be tested
n_ini_ply_drops = 10
# Minimum ply count for ply groups during ply drop layout generation
group_size_min = 8
# Desired ply count for ply groups during ply drop layout generation
group_size_max = 12
# Time limit to create a group ply-drop layout
time_limit_group_pdl = 1
# Time limit to create a ply-drop layout
time_limit_all_pdls = 10

## For step 3 of BELLA

# Branching limit for global pruning
global_node_limit = 100
# Branching limit for local pruning
local_node_limit = 100
# Branching limit for global pruning at the last level 
global_node_limit_final = 1
# Branching limit for local pruning at the last level 
local_node_limit_final = 100

## For step 4.1 of BELLA

# to save repair success rates
save_success_rate = True
# Thickness of the reference panels
n_plies_ref_panel = 1
# repair to improve the convergence towards the in-plane lamination parameter
# targets
repair_membrane_switch = True
# repair to improve the convergence towards the out-of-plane lamination
# parameter targets
repair_flexural_switch = True
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
global_node_limit2 = 5
# Branching limit for local pruning during ply drop layout optimisation
local_node_limit2 = global_node_limit2

## For step 4.3 of BELLA

# Branching limit for global pruning during ply drop layout optimisation
global_node_limit3 = 5
# Branching limit for local pruning during ply drop layout optimisation
local_node_limit3 = global_node_limit3

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
    time_limit_all_pdls=time_limit_all_pdls,
    save_buckling=True)

### Material properties -------------------------------------------------------

data_materials = pd.read_excel(filename_input, sheet_name='Materials',
                               header=None, index_col=0).T
E11 = data_materials["E11"].iloc[0]
E22 = data_materials["E22"].iloc[0]
nu12 = data_materials["nu12"].iloc[0]
G12 = data_materials["G12"].iloc[0]
density_area = data_materials["areal density"].iloc[0]
ply_t = data_materials["ply thickness"].iloc[0]
materials = Material(E11=E11, E22=E22, G12=G12, nu12=nu12,
                    density_area=density_area, ply_t=ply_t)

### Objective function parameters ---------------------------------------------

data_objective = pd.read_excel(filename_input, sheet_name='Objective function',
                               header=None, index_col=0).T
coeff_10 = data_objective["coeff_10"].iloc[0]
coeff_contig = data_objective["coeff_contig"].iloc[0]
coeff_diso = data_objective["coeff_diso"].iloc[0]
coeff_oopo = data_objective["coeff_oopo"].iloc[0]
coeff_spacing = data_objective["coeff_spacing"].iloc[0]

obj_func_param = ObjFunction(
    constraints=constraints,
    coeff_contig=coeff_contig,
    coeff_diso=coeff_diso,
    coeff_10=coeff_10,
    coeff_oopo=coeff_oopo,
    coeff_spacing=coeff_spacing)

### Multi-panel composite laminate layout -------------------------------------

data_panels = pd.read_excel(filename_input, sheet_name='Panels')

lampam_weightings_all = data_panels[[
    "lampam_weightings[1]", "lampam_weightings[2]", "lampam_weightings[3]",
    "lampam_weightings[4]", "lampam_weightings[5]", "lampam_weightings[6]",
    "lampam_weightings[7]", "lampam_weightings[8]", "lampam_weightings[9]",
    "lampam_weightings[10]", "lampam_weightings[11]", "lampam_weightings[12]"]]

lampam_targets_all = data_panels[[
    "lampam_target[1]", "lampam_target[2]", "lampam_target[3]",
    "lampam_target[4]", "lampam_target[5]", "lampam_target[6]",
    "lampam_target[7]", "lampam_target[8]", "lampam_target[9]",
    "lampam_target[10]", "lampam_target[11]", "lampam_target[12]"]]

panels = []
for ind_panel in range(data_panels.shape[0]):
    panels.append(Panel(
        ID=int(data_panels["Panel ID"].iloc[ind_panel]),
        lampam_target=np.array(lampam_targets_all.iloc[ind_panel], float),
        lampam_weightings=np.array(lampam_weightings_all.iloc[ind_panel], float),
        n_plies=int(data_panels["Number of plies"].iloc[ind_panel]),
        weighting=float(data_panels[
            "Weighting in MP objective funtion"].iloc[ind_panel]),
        neighbour_panels=np.array(data_panels[
            "Neighbour panel IDs"].iloc[ind_panel].split(" "), int),
        constraints=constraints,
        length_x=float(data_panels["Length_x"].iloc[ind_panel]),
        length_y=float(data_panels["Length_y"].iloc[ind_panel]),
        N_x=float(data_panels["N_x"].iloc[ind_panel]),
        N_y=float(data_panels["N_y"].iloc[ind_panel])))


multipanel = MultiPanel(panels)


### Optimiser Run -------------------------------------------------------------
result = BELLA_optimiser(multipanel, parameters, obj_func_param, constraints,
                         filename_result, materials)