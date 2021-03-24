# -*- coding: utf-8 -*-
"""
Class for reduced multi-panel structures

"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np

class ReducedMultiPanel():
    """
    Class for reduced multi-panel structures
    """
    def __init__(self, multipanel, constraints, n_plies_ref_panel=1):
        """Create object for storing multi-panel structures information

        - n_plies_in_panels:
            the ordered list with the panels ply counts that are used to build
            the stacking sequence table.
        - n_panels:
            number of panels with diffrent thicknesses in the stacking sequence
            table
        - n_plies_ref_panel:
            ply count of reference panel
        - ind_panels_guide:
            a mapping of the panels of the structure with the corrresponding
            stack in the stacking sequence table
        - ind_for_reduc:
            indices of n_panels panels with different thickness
        - ind_ref:
            index of reference panel for repair
        - boundaries:
            list of panels adjacency to consider expressed with the indices for
            the reduced number of panels
        - panel_weightings_ini_2_guide:
            initial panel weightings in the reduced multipanel objective
            function
        - panel_weightings_ini_2A_guide:
            initial panel weightings in the reduced multipanel in-plane
            objective function
        - panel_weightings_ini_2D_guide:
            initial panel weightings in the reduced multipanel out-of-plane
            objective function
        """

        (self.n_plies_in_panels, self.ind_for_reduc,
         self.ind_panels_guide) = np.unique(
             [panel.n_plies for panel in multipanel.panels],
             return_inverse=True, return_index=True)

        if constraints.sym:
            has_middle_ply = False
            for n_plies in self.n_plies_in_panels:
                if has_middle_ply and not n_plies % 2:
                    raise Exception("""
The panel ply counts are not compatible with guide-based blending (middle ply)""")
                if n_plies % 2:
                    has_middle_ply = True

        self.ind_for_reduc = self.ind_for_reduc.astype(int)

        self.n_panels = self.n_plies_in_panels.size
        self.ind_thick = self.n_panels - 1

#        print('n_plies_in_panels', self.n_plies_in_panels)
#        print('n_panels', self.n_panels)
#        print('ind_panels_guide', self.ind_panels_guide)
#        print('ind_for_reduc', self.ind_for_reduc)

        # reduced panel weightings in the reduced multipanel objective function
        self.panel_weightings_ini = np.zeros((self.n_panels,))
        for ind_panel in range(self.n_panels):
            self.panel_weightings_ini[self.ind_panels_guide[
                ind_panel]] += self.panel_weightings_ini[ind_panel]
#        print('self.panel_weightings_ini', self.panel_weightings_ini)

        if n_plies_ref_panel in self.n_plies_in_panels:
            ind_min = np.argmin(abs(self.n_plies_in_panels-n_plies_ref_panel))
            self.ind_ref = ind_min
        if n_plies_ref_panel <= self.n_plies_in_panels[0]:
            self.ind_ref = 0
        elif n_plies_ref_panel >= self.n_plies_in_panels[-1]:
            self.ind_ref = self.n_panels - 1
        else:
            index = 0
            while True:
                n_plies_1 = self.n_plies_in_panels[index]
                n_plies_2 = self.n_plies_in_panels[index + 1]
                if n_plies_ref_panel > n_plies_1 \
                and n_plies_ref_panel <= n_plies_2:
                    self.ind_ref = index + 1
                    break
                else:
                    index += 1

        self.n_plies_ref_panel = self.n_plies_in_panels[self.ind_ref]

        self.panels = [multipanel.panels[self.ind_for_reduc[ind_panel]] \
                       for ind_panel in range(self.n_panels)]

        self.n_panels_thin = self.ind_ref + 1
        self.n_panels_thick = self.n_panels - self.ind_ref

        self.middle_ply_indices = np.array(
            [self.panels[ind_panel].middle_ply_index \
             for ind_panel in range(self.n_panels)])

        self.check_boundary_weights(multipanel)

        self.calc_panel_weigthings(multipanel)

        self.calc_parameters(multipanel, constraints)

    def check_boundary_weights(self, multipanel):
        """
        check the boundary weights

        """
        self.boundaries = []
        for elem1, elem2 in multipanel.boundaries:
            elem1 = self.ind_panels_guide[elem1]
            elem2 = self.ind_panels_guide[elem2]
            loc_boundary = [elem1, elem2]
            loc_boundary.sort()
            if loc_boundary not in self.boundaries and elem1 != elem2:
                self.boundaries.append(loc_boundary)
        self.boundaries = np.array(self.boundaries)

        ## boundary weightings
        self.boundary_weights = dict()
        
        for panel1, panel2 in multipanel.boundaries:
            weight = multipanel.boundary_weights[(panel1, panel2)]
            panel1_in_strip = self.ind_panels_guide[panel1]
            panel2_in_strip = self.ind_panels_guide[panel2]
            panel1_in_strip, panel2_in_strip = sorted((panel1_in_strip, 
                                                       panel2_in_strip))
            if panel1_in_strip != panel2_in_strip:
            
                if (panel1_in_strip, panel2_in_strip) in self.boundary_weights:
                    self.boundary_weights[
                        (panel1_in_strip, panel2_in_strip)] += weight
                else:
                    self.boundary_weights[
                        (panel1_in_strip, panel2_in_strip)] = weight 
            

    def calc_panel_weigthings(self, multipanel):
        """
        returns the panel weightings for the objective function at each step
        of the thick-to-thin or thin-to-thick repair

        The function accounts only for the panel for which the ply drop will be
        determined.

        Moreover, the weightings of panel assumed to have the same lamination
        parameters are aggregated
        """
        self.actual_panel_weightings = np.zeros((self.n_panels,))
        for ind_panel, panel in enumerate(self.panels):
            self.actual_panel_weightings[
                self.ind_panels_guide[ind_panel]] += panel.weighting

        list_all = [[] for ind in range(self.n_panels)]

        for ind_panel in range(self.ind_ref):

    #        print('multipanel.reduced.actual_panel_weightings',
    #              multipanel.reduced.actual_panel_weightings)
            panel_weightings = np.copy(
                self.actual_panel_weightings)[:self.ind_ref]
            panel_weightings /= sum(panel_weightings)
    #        print('panel_weightings', panel_weightings)
            to_add, panel_weightings = (panel_weightings[:ind_panel],
                                        panel_weightings[ind_panel:])
    #        print('panel_weightings', panel_weightings, 'to_add', to_add)
            panel_weightings[0] += sum(to_add)
            list_all[ind_panel] = panel_weightings
    #        print('panel_weightings', panel_weightings)

        for ind_panel in range(self.ind_ref + 1, self.n_panels):

    #        print('self.actual_panel_weightings',
    #              self.actual_panel_weightings)
            panel_weightings = np.copy(
                self.actual_panel_weightings)[self.ind_ref + 1:] # YEEEES
            panel_weightings /= sum(panel_weightings)
    #        print('panel_weightings', panel_weightings)
            panel_weightings, to_add = (
                panel_weightings[:ind_panel - self.ind_ref],
                panel_weightings[ind_panel - self.ind_ref:])
    #        print('panel_weightings', panel_weightings, 'to_add', to_add)
            panel_weightings[-1] += sum(to_add)
            list_all[ind_panel] = panel_weightings
    #        print('panel_weightings', panel_weightings)

        self.panel_weightings = list_all


    def calc_parameters(self, multipanel, constraints):
        """
        caluclates values usefull during the beam search for thick-to-thin or
        thin-to-thick repair

        INPUTS

        - multipanel: multi-panel structure

        OUTPUTS

        - n_steps: number of ply drops to investigate
        - ind_panel_tab: index of panel where the next ply drop is located at each
            step of the beam search
        - n_plies_panel_after_tab: total number of plies in panel where the next
            ply drop is located at each step of the beam search
        - n_plies_after_tab: number of plies in panel after the next ply drop dealt
            with at each step of the beam search
        - n_plies_before_tab: number of plies in panel before the next ply drop
            dealt with at each step of the beam search
        - new_boundary_tab: indicate if the next ply drop to optimise is the first
            to be optimised in the panel, at each step of the beam search
        - n_panels_designed_tab: number of panels for which the ply drops have been
            designed so far
        """
        ### thick to thin
        n_steps = self.n_plies_ref_panel - self.n_plies_in_panels[0]

        if constraints.sym:
            n_steps //= 2

        if not n_steps:
            self.n_steps_thin = 0
            self.ind_panel_thin_tab = None
            self.n_plies_after_thin_tab = None
            self.n_plies_before_thin_tab = None
            self.new_boundary_thin_tab = None
            self.n_panels_designed_thin_tab = None

        else:

            ind_panel_tab = np.zeros((n_steps,), dtype='int16')
            n_plies_panel_after_tab = np.zeros((n_steps,), dtype='int16')
            n_plies_after_tab = np.zeros((n_steps,), dtype='int16')
            n_plies_before_tab = np.zeros((n_steps,), dtype='int16')
            n_panels_designed_tab = np.zeros((n_steps,), dtype='int16')
            new_boundary_tab = np.zeros((n_steps,), bool)

            ind_panel_tab[0] = self.ind_ref - 1
            n_plies_panel_after_tab[0] = self.n_plies_in_panels[
                ind_panel_tab[0]]
            if constraints.sym:
                n_plies_after_tab[0] = self.n_plies_ref_panel - 2
            else:
                n_plies_after_tab[0] = self.n_plies_ref_panel - 1
            n_plies_before_tab[0] = self.n_plies_ref_panel
            new_boundary_tab[0] = True

            for ind_step in range(1, n_steps):
                # add same values
                ind_panel_tab[ind_step] = ind_panel_tab[ind_step - 1]
                n_plies_after_tab[ind_step] = n_plies_after_tab[ind_step - 1]
                n_plies_before_tab[ind_step] = n_plies_before_tab[ind_step - 1]
                n_plies_panel_after_tab[ind_step] = n_plies_panel_after_tab[
                    ind_step - 1]

                if constraints.sym:
                    n_plies_after_tab[ind_step] -= 2
                    n_plies_before_tab[ind_step] -= 2
                else:
                    n_plies_after_tab[ind_step] -= 1
                    n_plies_before_tab[ind_step] -= 1
                if n_plies_after_tab[ind_step] < n_plies_panel_after_tab[ind_step]:
                    ind_panel_tab[ind_step] -= 1
                    n_plies_panel_after_tab[
                        ind_step] = self.n_plies_in_panels[
                            ind_panel_tab[ind_step]]
                    new_boundary_tab[ind_step] = True
                else:
                    new_boundary_tab[ind_step] = False

                n_panels_designed_tab[ind_step] = n_panels_designed_tab[ind_step - 1]
                if new_boundary_tab[ind_step]:
                    n_panels_designed_tab[ind_step] += 1

            self.n_steps_thin = n_steps
            self.ind_panel_thin_tab = ind_panel_tab
            self.n_plies_after_thin_tab = n_plies_after_tab
            self.n_plies_before_thin_tab = n_plies_before_tab
            self.new_boundary_thin_tab = new_boundary_tab
            self.n_panels_designed_thin_tab = n_panels_designed_tab

        ### thin to thick

        n_steps = self.n_plies_in_panels[-1] - self.n_plies_ref_panel

        if constraints.sym:
            n_steps //= 2

        if not n_steps:
            self.n_steps_thick = 0
            self.ind_panel_thick_tab = None
            self.n_plies_after_thick_tab = None
            self.n_plies_before_thick_tab = None
            self.new_boundary_thick_tab = None
            self.n_panels_designed_thick_tab = None
        else:
            ind_panel_tab = np.zeros((n_steps,), dtype='int16')
            n_plies_panel_after_tab = np.zeros((n_steps,), dtype='int16')
            n_plies_after_tab = np.zeros((n_steps,), dtype='int16')
            n_plies_before_tab = np.zeros((n_steps,), dtype='int16')
            n_panels_designed_tab = np.zeros((n_steps,), dtype='int16')
            new_boundary_tab = np.zeros((n_steps,), bool)

            ind_panel_tab[0] = self.ind_ref + 1
            n_plies_panel_after_tab[0] = self.n_plies_in_panels[
                ind_panel_tab[0]]
            if constraints.sym:
                n_plies_after_tab[0] = self.n_plies_ref_panel + 2
            else:
                n_plies_after_tab[0] = self.n_plies_ref_panel + 1
            n_plies_before_tab[0] = self.n_plies_ref_panel
            new_boundary_tab[0] = True

            for ind_step in range(1, n_steps):
                # add same values
                ind_panel_tab[ind_step] = ind_panel_tab[ind_step - 1]
                n_plies_after_tab[ind_step] = n_plies_after_tab[ind_step - 1]
                n_plies_before_tab[ind_step] = n_plies_before_tab[ind_step - 1]
                n_plies_panel_after_tab[ind_step] = n_plies_panel_after_tab[
                    ind_step - 1]

                if constraints.sym:
                    n_plies_after_tab[ind_step] += 2
                    n_plies_before_tab[ind_step] += 2
                else:
                    n_plies_after_tab[ind_step] += 1
                    n_plies_before_tab[ind_step] += 1
                if n_plies_after_tab[ind_step] > n_plies_panel_after_tab[ind_step]:
                    ind_panel_tab[ind_step] += 1
                    n_plies_panel_after_tab[
                        ind_step] = self.n_plies_in_panels[
                            ind_panel_tab[ind_step]]
                    new_boundary_tab[ind_step] = True
                else:
                    new_boundary_tab[ind_step] = False

                n_panels_designed_tab[ind_step] = n_panels_designed_tab[ind_step - 1]
                if new_boundary_tab[ind_step]:
                    n_panels_designed_tab[ind_step] += 1


            self.n_steps_thick = n_steps
            self.ind_panel_thick_tab = ind_panel_tab
            self.n_plies_after_thick_tab = n_plies_after_tab
            self.n_plies_before_thick_tab = n_plies_before_tab
            self.new_boundary_thick_tab = new_boundary_tab
            self.n_panels_designed_thick_tab = n_panels_designed_tab

    def __repr__(self):
        " Display object "
        return f"""
Reduced multipanel structure (blending strip):
    Number of panels: {self.n_panels}
    Number of plies per panel (ordered): {self.n_plies_in_panels}
    Number of plies in reference panel: {self.n_plies_ref_panel}
    Indices of panel middle plies: {self.middle_ply_indices}
    Reduced index of thick panel: {self.ind_thick}
    Reduced index of reference panel: {self.ind_ref}
    Panel boundaries in reduced indices: {self.boundaries}
    Mapping multipanel structure to blending strip: {self.ind_panels_guide}
    Mapping blending strip to multipanel structure: {self.ind_for_reduc}
"""

if __name__ == "__main__":

    import sys
    sys.path.append(r'C:\BELLA')
    from src.BELLA.panels import Panel
    from src.BELLA.constraints import Constraints
    from src.BELLA.multipanels import MultiPanel
    constraints = Constraints(sym=True)
    panel1 = Panel(1, constraints, neighbour_panels=[5], n_plies=12)
    panel2 = Panel(5, constraints, neighbour_panels=[1], n_plies=16)
    panel3 = Panel(2, constraints, neighbour_panels=[1, 5], n_plies=14)
    panel4 = Panel(6, constraints, neighbour_panels=[1], n_plies=10)
    multipanel = MultiPanel([panel1, panel2, panel3, panel4])
    print(multipanel)
    multipanel.from_mp_to_blending_strip(constraints, n_plies_ref_panel=12)
    print(multipanel.reduced)
#    print(multipanel.reduced.panels)
