# -*- coding: utf-8 -*-
"""
This script generates input ply drop layouts
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import pandas as pd

sys.path.append(r'C:\BELLA')
from src.BELLA.panels import Panel
from src.BELLA.multipanels import MultiPanel
from src.BELLA.parameters import Parameters
from src.BELLA.obj_function import ObjFunction
from src.BELLA.constraints import Constraints
from src.BELLA.divide_panels import divide_panels
from src.BELLA.pdl_ini import create_initial_pdls
from src.divers.excel import autofit_column_widths
from src.divers.excel import delete_file
from src.divers.excel import append_df_to_excel
from src.BELLA.save_set_up import save_constraints_BELLA
from src.BELLA.save_set_up import save_parameters_BELLA
from src.guidelines.ply_drop_spacing import calc_penalty_spacing
#from src.guidelines.one_stack import check_ply_drop_rules

# Number of initial ply drops to be tested
n_ini_ply_drops = 6
filename = 'pdls_ini_6_panels_5_boundaries.xlsx'
filename = 'pdls_ini_6_panels_9_boundaries.xlsx'
delete_file(filename)
#==============================================================================
# Targets and panel geometries
#==============================================================================
# panel number of plies
n_plies = [60, 56, 52, 48, 44, 40]
# number of panels
n_panels = len(n_plies)
# panel IDs
ID = list(range(1, n_panels + 1))
# panels adjacency
# panels adjacency
neighbour_panels = {1:[2],
                    2:[1, 3],
                    3:[2, 4],
                    4:[3, 5],
                    5:[4, 6],
                    6:[5]}
neighbour_panels = {1:[2, 3],
                    2:[1, 3, 4],
                    3:[2, 4, 1, 5],
                    4:[3, 5, 2, 6],
                    5:[3, 4, 6],
                    6:[5, 4]}
#==============================================================================
# Design guidelines
#==============================================================================
# symmetry
sym = True

# damage tolerance
dam_tol = True

# covering
covering = False

# ply drop spacing
pdl_spacing = True
# Minimum number of continuous plies required between two blocks of dropped
# plies
min_drop = 5

constraints = Constraints(
    sym=sym,
    dam_tol=dam_tol,
    covering=covering,
    min_drop=min_drop,
    pdl_spacing=pdl_spacing)



********* to modify *********

obj_func_param = ObjFunction(constraints)



#==============================================================================
# Optimiser Parameters
#==============================================================================
# Minimum ply count for ply groups during ply drop layout generation
group_size_min = 8
# Desired ply count for ply groups during ply drop layout generation
group_size_max = 12

# Coefficient for the ply drop spacing guideline penalty
coeff_spacing = 1

# Time limit to create a group ply-drop layout
time_limit_group_pdl = 1
# Time limit to create a ply-drop layout
time_limit_all_pdls = 100

# DO NOT DELETE
parameters = Parameters(
    constraints=constraints,
    group_size_min=group_size_min,
    group_size_max=group_size_max,
    n_ini_ply_drops=n_ini_ply_drops,
    coeff_spacing=coeff_spacing,
    time_limit_group_pdl=time_limit_group_pdl,
    time_limit_all_pdls=time_limit_all_pdls)

panels = []
for ind_panel in range(n_panels):
    panels.append(Panel(
        ID=ID[ind_panel],
        n_plies=n_plies[ind_panel],
        neighbour_panels=neighbour_panels[ID[ind_panel]],
        constraints=constraints))
#print(panels[0])

multipanel = MultiPanel(panels)
#print(multipanel)

#==============================================================================
# Ply drop layouts generartion
#==============================================================================
divide_panels(multipanel, parameters, constraints)
pdls = create_initial_pdls(multipanel, constraints, parameters, obj_func_param)

save_constraints_BELLA(filename, constraints)
save_parameters_BELLA(filename, parameters)

for ind_pdl, pdl in enumerate(pdls):
    table_pdl = pd.DataFrame()
    for ind_row, pdl_row in enumerate(pdl):
        for ind_elem, elem in enumerate(pdl_row):
            table_pdl.loc[ind_row, str(ind_elem)] = elem
    append_df_to_excel(
        filename, table_pdl, 'pdl' + str(ind_pdl+1), index=False, header=False)
autofit_column_widths(filename)

#==============================================================================
# Save ply drop penalties
#==============================================================================
table_penalties = pd.DataFrame()
for ind_pdl, pdl in enumerate(pdls):
    penalty_spacing = calc_penalty_spacing(
        pdl=pdl,
        multipanel=multipanel,
        constraints=constraints,
        on_blending_strip=True

        **are you sure?)

    table_penalties.loc[ind_pdl, 'penalty_spacing'] = penalty_spacing
    table_penalties.loc[ind_pdl, 'min_drop'] = constraints.min_drop

append_df_to_excel(
    filename, table_penalties, 'penalties', index=True, header=True)
autofit_column_widths(filename)
