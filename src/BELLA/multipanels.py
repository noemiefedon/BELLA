# -*- coding: utf-8 -*-
"""
Class for a multi-panel structure

"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np
sys.path.append(r'C:\BELLA')
from src.BELLA.parameters import Parameters
from src.BELLA.constraints import Constraints
from src.BELLA.panels import Panel
from src.BELLA.reduced_multipanels import ReducedMultiPanel

class MultiPanel():
    """
    Class for multi-panel structures
    """
    def __init__(self, panels, boundary_weights=None):
        """Create object for storing multi-panel structures information"""

        # list of panels (classes)
        self.panels = panels
        if not isinstance(panels, list):
            raise MultiPanelDefinitionError(
                'Attention, panels must be a list!')

        # total area of structure
        self.area = sum([el.area for el in panels])

        # total area for all patches
        self.area_patches = sum([el.area * el.n_plies for el in panels])

        # minimum ply count
        self.n_plies_min = min((el.n_plies for el in panels))

        # maximum ply count
        self.n_plies_max = max((el.n_plies for el in panels))
        self.is_thick_panels = [panel.n_plies == self.n_plies_max \
                                  for panel in self.panels]

        # number of panels
        self.n_panels = len(panels)

        # number of plies in the laminates
        self.n_plies_in_panels = np.array([self.panels[ind_panel].n_plies \
                            for ind_panel in range(self.n_panels)])

        self.has_a_middle_ply()
        self.identify_one_thickest_panel()

        self.calc_panel_boundary_dict(panels, boundary_weights)

    def should_you_use_BELLA(self):
        """ Tells the user when using LAYLA is better than employing BELLA

        Displays a message when BELLA is employed to design a composite
        laminate structure with one panel to indicate that LAYLA is better
        suited for the task than BELLA.

        Returns
        -------
        None.

        """

        if self.n_panels == 1:
            print("""
You are using BELLA to design a composite laminate structure with one panel.
LAYLA is better suited for this task than BELLA, please consider using LAYLA
instead of BELLA.""")

    def filter_target_lampams(self, constraints, obj_func_param):
        """
        filters applied to the lamination parameters to account for orthotropy
        requirements
        """
        for panel in self.panels:
            panel.filter_target_lampams(constraints, obj_func_param)

    def filter_lampam_weightings(self, constraints, obj_func_param):
        """
        filter of the lamination-parameter weightings in the panel
        objective function to account for the design guidelines

        lampam_weightings_3: for blending steps 3 (contain penalty for
            out-of-plane orthotropy and may contain penalty for balance)
        lampam_weightings: for all other blending steps (contain penalty for
            out-of-plane orthotropy and does not contain penalty for balance)
        """
        for panel in self.panels:
            panel.filter_lampam_weightings(constraints, obj_func_param)

    def from_mp_to_blending_strip(self, constraints, n_plies_ref_panel=1):
        """
        performs the blending step 2: maps the multi-panel structure to a
        blending strip, i.e. a series a panels in a row
        """
        self.reduced = ReducedMultiPanel(self, constraints, n_plies_ref_panel)


    def calc_panel_boundary_dict(self, panels, boundary_weights):
        """
        checks that all panels have a different ID
        collates all the panel boundaries in self.boundaries
        checks that all panels are connected
        """
        ## checks that all panels have a different ID
        self.dict_ID_to_indices = dict()
        for ind_panel, panel in enumerate(panels):
            panel.ID_code = ind_panel
            self.dict_ID_to_indices[panel.ID] = ind_panel
        if len(self.dict_ID_to_indices) != self.n_panels:
            raise MultiPanelDefinitionError("""
Several panels with the same index!""")
#        print('dict_ID_to_indices', self.dict_ID_to_indices)

        ## create the dictionary of panel boundaries
        self.boundaries = []
        for ind_panel, panel in enumerate(panels):
            neighbours = [self.dict_ID_to_indices[neighbour] \
                          for neighbour in panel.neighbour_panels]
            for elem in neighbours:
                self.boundaries.append(np.sort([ind_panel, elem]))
                self.boundaries.append(np.flip(np.sort([ind_panel, elem])))
        if len(self.boundaries) == 0:
            self.boundaries = np.array((), int).reshape((0,2))
        else:
            self.boundaries = np.unique(self.boundaries, axis=0)
#        print('boundaries', self.boundaries)

        ## checks that all panels are connected
        visited_nodes = []
        set_avail_nodes = set([0])
        while len(set_avail_nodes) != 0 and len(visited_nodes) < self.n_panels:
            current_node = set_avail_nodes.pop()
            visited_nodes.append(current_node)
            for elem in self.boundaries:
                if elem[0] == current_node and elem[1] not in visited_nodes\
                and elem[1] not in set_avail_nodes:
                    set_avail_nodes.add(elem[1])
#        print('visited_nodes', visited_nodes)
        if not len(visited_nodes) == self.n_panels:
            raise MultiPanelDefinitionError("""
The panels of the multipanel-component are not all connected!""")

        if len(self.boundaries) == 0:
            self.boundaries = np.array((), int).reshape((0,2))
        else:
            self.boundaries = np.unique(
                np.array([np.sort(elem) for elem in self.boundaries]), axis=0)
#        print('boundaries', self.boundaries)

        ## dictionary with panel Ids
        self.boundaries_in_IDs = np.empty((self.boundaries.shape[0], 2), int)
        for ind_row, (first, second) in enumerate(self.boundaries):
            self.boundaries_in_IDs[ind_row, 0] = self.panels[first].ID
            self.boundaries_in_IDs[ind_row, 1] = self.panels[second].ID


        ## reorganise the boundary weightings
        self.boundary_weights_in_IDs = dict()
        self.boundary_weights = dict()
        if boundary_weights is not None:

            for weight in boundary_weights.values():
                if weight < 0:
                    raise Exception(
                        'The boundary weightings should be positive.')

            if len(boundary_weights) < self.boundaries.shape[0]:
                print(len(boundary_weights), self.boundaries)
                raise Exception(
                    'Insufficient number of boundary weightings.')

            for ind_panel1, ind_panel2 in self.boundaries_in_IDs:
                ind_panel1_mod = self.dict_ID_to_indices[ind_panel1]
                ind_panel2_mod = self.dict_ID_to_indices[ind_panel2]
                ind_panel1, ind_panel2 = sorted((ind_panel1, ind_panel2))
                ind_panel1_mod, ind_panel2_mod = sorted((ind_panel1_mod,
                                                         ind_panel2_mod))
                weight = boundary_weights.get((ind_panel1, ind_panel2), None)
                if weight:
                    self.boundary_weights_in_IDs[
                        (ind_panel1, ind_panel2)] = weight
                    self.boundary_weights[
                        (ind_panel1_mod, ind_panel2_mod)] = weight
                else:
                    weight = boundary_weights.get(
                        (ind_panel2, ind_panel1), None)
                    if not weight:
                        raise Exception('Missing boundary weightings.')
                    self.boundary_weights_in_IDs[
                        (ind_panel2, ind_panel1)] = weight
                    self.boundary_weights[
                        (ind_panel2_mod, ind_panel1_mod)] = weight

        else:  # all boundary weightings set to one
            for ind_panel1, ind_panel2 in self.boundaries_in_IDs:
                ind_panel1_mod = self.dict_ID_to_indices[ind_panel1]
                ind_panel2_mod = self.dict_ID_to_indices[ind_panel2]
                ind_panel1, ind_panel2 = sorted((ind_panel1, ind_panel2))
                ind_panel1_mod, ind_panel2_mod = sorted((ind_panel1_mod,
                                                         ind_panel2_mod))
                self.boundary_weights_in_IDs[(ind_panel1, ind_panel2)] = 1
                self.boundary_weights[(ind_panel1_mod, ind_panel2_mod)] = 1

        return 0

    def has_a_middle_ply(self):
        """
        returns:
            - middle_ply_indices: the locations of middle plies per panel
            (0 if no middle ply)
            - has_middle_ply: True if one panel at least has a middle ply
            - thick_panel_has_middle_ply: True if thickest panel has a middle
            ply
        """
        # locations of middle plies per panel (0 if no middle ply)
        self.middle_ply_indices = np.array(
            [self.panels[ind_panel].middle_ply_index \
             for ind_panel in range(self.n_panels)])
        self.has_middle_ply = bool(sum(self.middle_ply_indices))

        if self.has_middle_ply and self.n_plies_max % 2:
            self.thick_panel_has_middle_ply = True
        else:
            self.thick_panel_has_middle_ply = False


    def calc_ply_drops(self, inner_step):
        """
        returns a vector of the number of ply drops at each panel boundary of
        the blending strip for the inner_step-eme group of plies
        """
        n_ply_drops = np.zeros((self.reduced.n_panels,), dtype='int16')
        for index, panel in enumerate(self.reduced.panels):
            n_ply_drops[index] = self.reduced.n_plies_per_group[inner_step] \
                        - panel.n_plies_per_group[inner_step]
        return n_ply_drops

    def calc_weight(self, density_area):
        """
        returns the weight of the multipanel structure
        """
        return density_area*sum([panel.area*panel.n_plies \
                                      for panel in self.panels])

    def calc_weight_per_panel(self, density_area):
        """
        returns the weight of the multipanel structure per panel
        """
        self.weight_ref_per_panel = density_area * \
        np.array([panel.area*panel.n_plies for panel in self.panels])

    def calc_weight_from_sst(self, sst, density_area):
        """
        returns the weight of the multipanel structure from a stacking sequence
        table
        """
        return density_area*sum([panel.area * sum(sst[ind_panel] != -1) \
                                 for ind_panel,
                                 panel in enumerate(self.panels)])


    def identify_neighbour_panels(self):
        """
        returns the indices of the neighbouring panels for each panel
        """
        liste = []
        for ind_panel in range(self.n_panels):
            liste.append([])
        for boundary in self.boundaries:
            liste[boundary[0]].append(boundary[1])
            liste[boundary[1]].append(boundary[0])
        return liste


    def identify_one_thickest_panel(self):
        """
        returns the index of one of the thickest panels
        """
        for ind_panel, panel in enumerate(self.panels):
            if panel.n_plies == self.n_plies_max:
                self.ind_thick = ind_panel
                return 0
        raise Exception("""
The maximum number of plies should be the ply count of a panel""")


    def identify_thickest_panels(self, sym=False):
        """
        returns the index of all of the thickest panels
        """
        liste = []
        if sym and self.n_plies_max % 2 == 1: # midlle ply in thickest panels
            for ind_panel, panel in enumerate(self.panels):
                if panel.n_plies == self.n_plies_max \
                or panel.n_plies == self.n_plies_max - 1:
                    liste.append(ind_panel)
        else:
            for ind_panel, panel in enumerate(self.panels):
                if panel.n_plies == self.n_plies_max:
                    liste.append(ind_panel)
        if liste:
            return liste
        raise Exception("""
The maximum number of plies should be the ply count of a panel""")


    def __repr__(self):
        " Display object "

        to_add = ''
        # number of groups
        if hasattr(self, 'n_groups'):
            to_add = to_add + 'Number of groups : ' +  str(self.n_groups) \
            + '\n'
        # number of plies per group for thickest laminates
        if hasattr(self, 'n_plies_per_group'):
            to_add = to_add + 'Max number of plies per group : ' \
                  +  str(self.n_plies_per_group) + '\n'
        # position of the group first plies for thickest laminates
        if hasattr(self, 'n_first_plies'):
            to_add = to_add + 'Position first plies : ' \
                  +  str(self.n_first_plies) + '\n'

        return f"""
Number of panels : {self.n_panels}
Maximum number of plies in a panel: {self.n_plies_max}
Index of one of the thickest panels: {self.ind_thick}
Area : {self.area}
Area for all patches: {self.area_patches}
Panel boundary matrix : {self.boundaries_in_IDs}
""" + to_add


class MultiPanelDefinitionError(Exception):
    " Errors during the definition of a multi-panel structure"

if __name__ == "__main__":
    print('*** Test for the class MultiPanel ***\n')
    constraints = Constraints(
        sym=True,
        dam_tol=False,
        covering=False,
        pdl_spacing=True,
        min_drop=2)
    parameters = Parameters(constraints=constraints, n_plies_ref_panel=48)
    n_plies_target1 = 48
    n_plies_target2 = 46
    n_plies_target3 = 40
    n_plies_target4 = 40
    panel1 = Panel(ID=1,
                    n_plies=n_plies_target1,
                    constraints=constraints,
                    neighbour_panels=[2])
    panel2 = Panel(ID=2,
                    n_plies=n_plies_target2,
                    constraints=constraints,
                    neighbour_panels=[1, 3])
    panel3 = Panel(ID=3,
                    n_plies=n_plies_target3,
                    constraints=constraints,
                    neighbour_panels=[2, 4])
    panel4 = Panel(ID=4,
                    n_plies=n_plies_target4,
                    constraints=constraints,
                    neighbour_panels=[3])
    multipanel = MultiPanel([panel1, panel2, panel3, panel4])
    print(multipanel)

    from src.BELLA.divide_panels import divide_panels
    divide_panels(multipanel, parameters, constraints)

    print('multipanel.reduced.n_plies_in_panels', multipanel.reduced.n_plies_in_panels)
    print('multipanel.calc_ply_drops(0)', multipanel.calc_ply_drops(0))
    print('multipanel.reduced.n_plies_per_group', multipanel.reduced.n_plies_per_group)
    print('multipanel.reduced.middle_ply_indices', multipanel.reduced.middle_ply_indices)