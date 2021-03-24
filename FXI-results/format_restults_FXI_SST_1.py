# -*- coding: utf-8 -*-
"""
This script formats the results of blending optimisations the optimiser of
Francois-Xavier Irisarri
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
from src.BELLA.results import BELLA_Results
from src.BELLA.results import BELLA_ResultsOnePdl
from src.BELLA.format_pdl import convert_sst_to_ss
from src.guidelines.ipo_oopo import calc_penalty_oopo_ss
from src.guidelines.ipo_oopo import calc_penalty_ipo
from src.guidelines.contiguity import calc_penalty_contig_mp
from src.guidelines.disorientation import calc_number_violations_diso_mp
from src.guidelines.ten_percent_rule import calc_penalty_10_pc
from src.guidelines.ten_percent_rule import calc_ply_counts
from src.guidelines.ten_percent_rule import calc_penalty_10_ss
from src.guidelines.ply_drop_spacing  import calc_penalty_spacing
from src.BELLA.save_set_up import save_constraints_BELLA
from src.BELLA.save_set_up import save_multipanel, save_objective_function_BELLA
from src.BELLA.save_set_up import save_materials
from src.BELLA.save_result import save_result_BELLAs
from src.divers.excel import delete_file, autofit_column_widths

filename = 'SST-40-60.xlsx'
# filename = 'SST-80-120.xlsx'
# filename = 'SST-120-180.xlsx'
filename_input = '/BELLA/input-files/input_file_' + filename
filename_FXI = '/BELLA/FXI-results/result_FXI_' + filename
filename_res = 'results_FXI_' + filename

# check for authorisation before overwriting
delete_file(filename_res)

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
        constraints=constraints))

multipanel = MultiPanel(panels)
multipanel.filter_lampam_weightings(constraints, obj_func_param)

### Organise the data structures of the results of FXI ------------------------


sst_all = pd.read_excel(filename_FXI, sheet_name='Best result FXI').fillna(-1)
sst_all = np.array(sst_all, int).T

dic_n_plies_sst = {}
for line in sst_all:
    n_plies = line[line != -1].size
    dic_n_plies_sst[n_plies] = line

sst = []

for panel in multipanel.panels:
    sst.append(dic_n_plies_sst[panel.n_plies])
sst = np.array(sst)

# remove unecessary -1
for ind in range(sst.shape[1])[::-1]:
    if (sst[:, ind] == -1).all():
        sst = np.delete(sst, np.s_[ind], axis=1)

ss = convert_sst_to_ss(sst)

# lamination parameters
lampam = np.array([calc_lampam(ss[ind_panel]) \
                   for ind_panel in range(multipanel.n_panels)])

# disorientaion - penalty used in blending steps 4.2 and 4.3
n_diso = calc_number_violations_diso_mp(ss, constraints)
if constraints.diso and n_diso.any():
    penalty_diso = n_diso
else:
    penalty_diso = np.zeros((multipanel.n_panels,))

# contiguity - penalty used in blending steps 4.2 and 4.3
n_contig = calc_penalty_contig_mp(ss, constraints)
if constraints.contig and n_contig.any():
    penalty_contig = n_contig
else:
    penalty_contig = np.zeros((multipanel.n_panels,))

# 10% rule - no penalty used in blending steps 4.2 and 4.3
if constraints.rule_10_percent and constraints.rule_10_Abdalla:
    penalty_10 = calc_penalty_10_ss(ss, constraints, lampam, mp=True)
else:
    penalty_10 = calc_penalty_10_pc(
        calc_ply_counts(multipanel, ss, constraints), constraints)

# balance
penalty_bal_ipo = calc_penalty_ipo(lampam)

# out-of-plane orthotropy
penalty_oopo = calc_penalty_oopo_ss(lampam, constraints=constraints)

# penalty_spacing
penalty_spacing = calc_penalty_spacing(
    pdl=sst,
    multipanel=multipanel,
    constraints=constraints,
    on_blending_strip=False)

results = BELLA_Results(constraints, multipanel)
results_one_pdl = BELLA_ResultsOnePdl()

results_one_pdl.ss = ss
results_one_pdl.lampam = lampam
results_one_pdl.penalty_diso = penalty_diso
results_one_pdl.penalty_contig = penalty_contig
results_one_pdl.penalty_10 = penalty_10
results_one_pdl.penalty_bal_ipo = penalty_bal_ipo
results_one_pdl.penalty_oopo = penalty_oopo
results_one_pdl.penalty_spacing = penalty_spacing
results_one_pdl.n_diso = n_diso
results_one_pdl.n_contig = n_contig
results_one_pdl.sst = sst

results.update(0, results_one_pdl)
results.lampam = lampam
results.sst = sst
results.ss = ss

### Save data -----------------------------------------------------------------
save_constraints_BELLA(filename_res, constraints)
save_materials(filename_res, materials)
save_objective_function_BELLA(filename_res, obj_func_param)
save_multipanel(filename_res, multipanel, obj_func_param, materials)
save_result_BELLAs(filename_res, multipanel, constraints, None,
                   obj_func_param, None, results, materials, only_best=True)
autofit_column_widths(filename_res)