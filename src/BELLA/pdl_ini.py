# -*- coding: utf-8 -*-
"""
Create ply drop layout

All stacking sequences are subsets from the stacking sequence of the thickest
panels.

- read_pdls_excel
    reads ply drop layout inputs in Excel file

- check_pdls_ini
    checks that the initial ply drop layouts have the correct numbers of panels
    and numbers of plies per panel

- create_initial_pdls
    creates ply drop layouts for a blended structure

- create_initial_pdl
    creates a ply drop layout for a blended structure

"""
import sys
import time
import pandas as pd
import numpy as np
import numpy.matlib

sys.path.append(r'C:\BELLA')
from src.BELLA.panels import Panel
from src.BELLA.multipanels import MultiPanel
from src.BELLA.parameters import Parameters
from src.BELLA.constraints import Constraints
from src.BELLA.obj_function import ObjFunction
from src.BELLA.divide_panels import divide_panels
from src.BELLA.pdl_group import randomly_pdl_guide
from src.BELLA.pdl_tools import global_pdl_from_local_pdl
from src.guidelines.ply_drop_spacing import calc_penalty_spacing
from src.divers.pretty_print import print_list_ss
#from src.guidelines.ply_drop_spacing  import indic_violation_ply_drop_spacing
#from src.guidelines.one_stack import check_ply_drop_rules

def read_pdls_excel(filename):
    """
    reads ply drop layout inputs in Excel file
    """
    pdls = []
    ind_pdl = 1
    while True:
        try:
            pdl = pd.read_excel(
                filename, sheet_name='pdl' + str(ind_pdl), dtype='int',
                header=None, )
        except:
            break
        finally:
            ind_pdl += 1
            pdls.append(np.array(pdl))
    ind_pdl -= 2
    print('    ' + str(ind_pdl) + ' intial ply-drop layouts read.')
    return pdls


def check_pdls_ini(multipanel, pdls_ini):
    """
    checks that the initial ply drop layouts have the correct numbers of panels
    and numbers of plies per panel
    """
    for pdl in pdls_ini:
        if pdl.shape[0] != multipanel.reduced.n_panels:
            print(pdl.shape[0], multipanel.reduced.n_panels)
            raise Exception("""
Initial ply drop layout with wrong number of panels""")
        for ind_panel, pdl_panel in enumerate(pdl):
            if pdl_panel[pdl_panel != -1].size \
            != multipanel.reduced.panels[ind_panel].n_plies:
                raise Exception("""
Initial ply drop layout with wrong number of plies per panel""")
    print('    Sizes of the intial ply-drop layouts checked.')
    return None


def create_initial_pdls(multipanel, constraints, parameters, obj_func_param):
    """
    creates ply drop layouts for a blended structure with no duplicates
    """
    t_ini= time.time()

    pdls_perfect = np.zeros((
        parameters.n_ini_ply_drops,
        multipanel.reduced.n_panels,
        multipanel.n_plies_max), int)
    pdls_imperfect = np.zeros((
        parameters.n_ini_ply_drops,
        multipanel.reduced.n_panels,
        multipanel.n_plies_max), int)

    p_pdls_imperfect = np.zeros((parameters.n_ini_ply_drops,), dtype=float)

    ind_imperfect = 0
    ind_perfect = 0
    t_ini = time.time()
    elapsed_time = 0


    while ind_perfect < parameters.n_ini_ply_drops\
    and elapsed_time < parameters.time_limit_all_pdls:

#        print('ind_perfect', ind_perfect)
#        print('ind_imperfect', ind_imperfect)

        new_pdl = create_initial_pdl(multipanel, constraints, parameters,
                                     obj_func_param)

#        print('new_pdl')
#        print_list_ss(new_pdl)

#        if constraints.sym and multipanel.has_middle_ply:
#            pdl_after=np.flip(new_pdl[:, 1:], axis=1)
#        elif constraints.sym:
#            pdl_after=np.flip(new_pdl[:, :], axis=1)
#        else:
#            pdl_after=None

        new_penalty_spacing = calc_penalty_spacing(
            pdl=new_pdl,
            pdl_after=None,
            multipanel=multipanel,
            obj_func_param=obj_func_param,
            constraints=constraints,
            on_blending_strip=True)

        elapsed_time = time.time() - t_ini

        # Store the new pdl if it is perfect (no violation of manufacturing
        # constraint) or if it is among the parameters.n_ini_ply_drops best
        # unmanufacturable solutions found so far
        if new_penalty_spacing == 0:

            # To remove duplicates
            is_double = False
            for ind in range(ind_perfect):
                if (new_pdl - pdls_perfect[ind] == 0).all():
                    is_double = True
                    break
            if is_double:
                continue
            pdls_perfect[ind_perfect] = new_pdl
            ind_perfect += 1

        else:
            # To only keep the imperfect pdl with the smallest penalties
            if ind_imperfect >= parameters.n_ini_ply_drops:
                if new_penalty_spacing < max(p_pdls_imperfect):
                    # To remove duplicates
                    is_double = False
                    for ind in range(ind_imperfect):
                        if np.allclose(new_pdl, pdls_imperfect[ind]):
                            is_double = True
                            break
                    if is_double:
                        continue
                    #print('is_double', is_double)
                    indexx = np.argmin(p_pdls_imperfect)
                    pdls_imperfect[indexx] = new_pdl
                    p_pdls_imperfect[indexx] = new_penalty_spacing
            else:
                # To remove duplicates
                is_double = False
                for ind in range(ind_imperfect):
                    if np.allclose(new_pdl, pdls_imperfect[ind]):
                        is_double = True
                        break
                if is_double:
                    continue
                #print('is_double', is_double)
#                print(new_pdl)
#                print(ind_imperfect)
#                print(pdls_imperfect.shape)
                pdls_imperfect[ind_imperfect] = new_pdl
                p_pdls_imperfect[ind_imperfect] = new_penalty_spacing
                ind_imperfect += 1

    # if the time limit is reached
    if elapsed_time >= parameters.time_limit_all_pdls:

        pdls_imperfect = pdls_imperfect[:ind_imperfect]
        p_pdls_imperfect = p_pdls_imperfect[:ind_imperfect]
#        print('pdls_perfect', pdls_perfect)
#        print('pdls_imperfect', pdls_imperfect)

        if not ind_imperfect + ind_perfect:
            raise Exception("""
No conform ply drop layout can be generated.
Too many ply drops between two adjacent panels.""")

        # add the non-manufacturable ply drop layouts
        for ind in range(ind_perfect, parameters.n_ini_ply_drops):
            indexx = np.argmin(p_pdls_imperfect)
            pdls_perfect[ind] = pdls_imperfect[indexx]
            p_pdls_imperfect[indexx] = 10e6

        print('    ' + str(ind_perfect) \
          + ' feasible intial ply-drop layouts generated.')
        print('    ' + str(parameters.n_ini_ply_drops - ind_perfect) \
          + ' infeasible intial ply-drop layouts generated.')

        return pdls_perfect

    # if enough manufacturable ply drop layouts have been found
    print('    ' + str(parameters.n_ini_ply_drops) \
          + ' feasible intial ply-drop layouts generated.')

    return pdls_perfect

def create_initial_pdl(multipanel, constraints, parameters, obj_func_param):
    """
    creates a ply drop layout for a blended structure

    - obj_func_param: objective function parameters
    """
    if constraints.sym:
        pdl_before_cummul = [None]*(multipanel.reduced.n_groups + 1)

        for index_pdl in range(len(pdl_before_cummul)):
            pdl_before_cummul[index_pdl] = None

        # plies for covering rule (including damage tolerance)
        if constraints.n_covering == 1:
            pdl_before_cummul[0] = np.matlib.repmat(
                np.array([0], dtype=int), multipanel.reduced.n_panels, 1)
        elif constraints.n_covering == 2:
            pdl_before_cummul[0] = np.matlib.repmat(
                np.array([0, 1], dtype=int), multipanel.reduced.n_panels, 1)

        for inner_step in range(multipanel.reduced.n_groups):

            last_group = bool(inner_step == multipanel.reduced.n_groups - 1)

            covering_top = False

            # create group ply drop layouts
            n_ply_drops = multipanel.calc_ply_drops(inner_step)
            my_pdl = randomly_pdl_guide(
                multipanel=multipanel,
                boundaries=multipanel.reduced.boundaries,
                has_middle_ply=multipanel.has_middle_ply,
                middle_ply_indices=multipanel.reduced.middle_ply_indices,
                n_ply_drops=n_ply_drops,
                n_max=multipanel.reduced.n_plies_per_group[inner_step],
                parameters=parameters,
                obj_func_param=obj_func_param,
                constraints=constraints,
                pdl_before=pdl_before_cummul[inner_step],
                last_group=last_group,
                covering_top=covering_top)
            pdl_before_cummul[inner_step + 1] = my_pdl[0]

        return global_pdl_from_local_pdl(
            multipanel, constraints.sym, pdl_before_cummul)

    pdl_before_cummul = [None]*(multipanel.reduced.n_groups + 2)
    pdl_after_cummul = [None]*(multipanel.reduced.n_groups + 2)
    for index_pdl in range(len(pdl_before_cummul)):
        pdl_before_cummul[index_pdl] = None
        pdl_after_cummul[index_pdl] = None

    if constraints.n_covering == 1:
        pdl_before_cummul[0] = np.matlib.repmat(
            np.array([0], dtype=int), multipanel.reduced.n_panels, 1)
        pdl_after_cummul[1] = np.matlib.repmat(
            np.array([0], dtype=int), multipanel.reduced.n_panels, 1)
    elif constraints.n_covering == 2:
        pdl_before_cummul[0] = np.matlib.repmat(
            np.array([0, 1], dtype=int), multipanel.reduced.n_panels, 1)
        pdl_after_cummul[1] = np.matlib.repmat(
            np.array([0, 1], dtype=int), multipanel.reduced.n_panels, 1)

    for inner_step in range(multipanel.reduced.n_groups):
        last_group = bool(inner_step == multipanel.reduced.n_groups - 1)

        if inner_step % 2 == 0:
            pdl_before = pdl_before_cummul[inner_step]
            pdl_after = None
        else:
            pdl_before = None
            pdl_after = pdl_after_cummul[inner_step]

        covering_top = False
        covering_bottom = False

        # create group ply drop layouts
        n_ply_drops = multipanel.calc_ply_drops(inner_step)
        my_pdl = randomly_pdl_guide(
            multipanel=multipanel,
            boundaries=multipanel.reduced.boundaries,
            n_ply_drops=n_ply_drops,
            n_max=multipanel.reduced.n_plies_per_group[inner_step],
            pdl_before=pdl_before,
            pdl_after=pdl_after,
            last_group=last_group,
            parameters=parameters,
            obj_func_param=obj_func_param,
            constraints=constraints,
            covering_top=covering_top,
            covering_bottom=covering_bottom)
#        print(my_pdl)
        if inner_step % 2 == 0:
            pdl_before_cummul[inner_step + 2] = my_pdl[0]
        else:
            pdl_after_cummul[inner_step + 2] = my_pdl[0]

    return global_pdl_from_local_pdl(
        multipanel, constraints.sym, pdl_before_cummul, pdl_after_cummul)


if __name__ == "__main__":
    print('\n*** Test for the function create_initial_pdl ***')
    constraints = Constraints(
        sym=True,
        dam_tol=False,
        covering=False,
        pdl_spacing=True,
        min_drop=2)
    parameters = Parameters(constraints=constraints)
    obj_func_param = ObjFunction(constraints)
    n_plies_target1 = 21
    n_plies_target2 = 19
    n_plies_target3 = 18
    n_plies_target4 = 16
    panel_1 = Panel(ID=1,
                    n_plies=n_plies_target1,
                    constraints=constraints,
                    neighbour_panels=[2])
    panel_2 = Panel(ID=2,
                    n_plies=n_plies_target2,
                    constraints=constraints,
                    neighbour_panels=[1, 3])
    panel_3 = Panel(ID=3,
                    n_plies=n_plies_target3,
                    constraints=constraints,
                    neighbour_panels=[2, 4])
    panel_4 = Panel(ID=4,
                    n_plies=n_plies_target4,
                    constraints=constraints,
                    neighbour_panels=[3])
    multipanel = MultiPanel(panels=[panel_1, panel_2, panel_3, panel_4])
    multipanel.from_mp_to_blending_strip(constraints)
    divide_panels(multipanel, parameters, constraints)
    pdl_ini = create_initial_pdl(
        multipanel,
        constraints,
        parameters,
        obj_func_param)
    print(pdl_ini)

