# -*- coding: utf-8 -*-
"""
Class for panels

"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np
sys.path.append(r'C:\BELLA')
from src.BELLA.parameters import Parameters
from src.BELLA.constraints import Constraints

class Panel():
    """
    Class for panels
    """
    def __init__(self,
                 ID,
                 constraints,
                 n_plies,
                 lampam_target=np.array([]),
                 lampam_weightings=np.ones((12,), float),
                 area=0,
                 length_x=0,
                 length_y=0,
                 N_x=0,
                 N_y=0,
                 weighting=1,
                 neighbour_panels=[]):
        """
        Create object for storing information concerning a panel:
        - n_plies: target number of plies
        - lampam_target: lamination-parameter targets
        - lampam_weightings: lamination-parameter weightings in panel objective
        function
        - weighting: panel weighting in the multi-panel objecive function
        - area: panel area
        - length_x: panel length (x-direction)
        - length_y: panel width (y-direction)
        - N_x: loading intensity in the x-direction
        - N_y: loading intensity in the y-direction
        - neighbour_panels: list of the neighbour panels' indices
        - constraints: design and manufacturing constraints
        - parameters: parameters of the optimiser
        """
        # set the number of plies and the number of a potential middle ply
        # for symmetric laminates
        self.set_n_plies(n_plies, constraints)

        # lamination-parameter targets
        self.lampam_target = lampam_target

        # panel ID
        self.ID = ID

        # list of the neighbour panels' indices
        self.neighbour_panels = neighbour_panels

        # panel area
        self.area = area
        self.length_x = length_x
        self.length_y = length_y
        if not isinstance(area, (float, int)):
            raise PanelDefinitionError("""
Attention, the panel area must be a number!
""")
        if area < 0:
            raise PanelDefinitionError("""
Attention, the panel area must be positive!
""")
        if not isinstance(length_x, (float, int)):
            raise PanelDefinitionError("""
Attention, the panel length (length_x) must be a number!
""")
        if length_x < 0:
            raise PanelDefinitionError("""
Attention, the panel length (length_x) must be positive!
""")
        if not isinstance(length_y, (float, int)):
            raise PanelDefinitionError("""
Attention, the panel length (length_y) must be a number!
""")
        if length_y < 0:
            raise PanelDefinitionError("""
Attention, the panel length (length_y) must be positive!
""")
        if length_x != 0 and length_y != 0:
            self.area = length_x * length_y

        # panel loading conditions
        self.N_x = N_x
        self.N_y = N_y

        # panel weighting in the multi-panel objective function
        self.weighting = weighting # initial value
        if not isinstance(weighting, (float, int)):
            raise PanelDefinitionError("""
The panel weighting in the multi-panel objective function must a number!
""")
        if weighting < 0:
            raise PanelDefinitionError("""
The panel weighting in the multi-panel objective function must be positive!
""")

        # lamination-parameter weightings in panel objective function
        if not(isinstance(lampam_weightings, np.ndarray)) \
        and lampam_weightings.size == 12 \
        and lampam_weightings.dtype == float:
            raise PanelDefinitionError("""
Attention, lampam_weightings must be a vector with 12 float components!
""")
        if [False for elem in lampam_weightings if elem < 0]:
                    raise PanelDefinitionError("""
Attention, the elements of lampam_weightings must be positive!
""")
        self.lampam_weightings_ini = lampam_weightings

    def filter_target_lampams(self, constraints, obj_func_param):
        """
        filters applied to the lamination parameters to account for orthotropy
        requirements
        """
        # If symmetry is desired, the corresponding target amination parameters
        # must be set to 0
        if constraints.sym:
            self.lampam_target[4:8] = 0
            self.lampam_target[4:8] = 0
        # If the balance rule is desired, the corresponding target
        # lamination parameters must be set to 0
        if constraints.bal:
            self.lampam_target[2] = 0
            self.lampam_target[3] = 0
        # If the out-of-plane orthotropy is desired, the corresponding target
        # lamination parameters must be set to 0
        if constraints.oopo:
            self.lampam_target[10] = 0
            self.lampam_target[11] = 0

    def filter_lampam_weightings(self, constraints, obj_func_param):
        """
        filter of the lamination-parameter weighting in the panel
        objective function to account for the design guidelines

#        lampam_weightings_3: for blending steps 3 (contain penalty for
#            out-of-plane orthotropy and may contain penalty for balance)
        lampam_weightings: for all other blending steps (contain penalty for
            out-of-plane orthotropy and does not contain penalty for balance)
        """
        lampam_weightings = np.copy(self.lampam_weightings_ini)

        # filter for zero lamination parameters factors
        if constraints.sym:
            lampam_weightings[4:8] = 0
        if set(constraints.set_of_angles) == set([0, 45, -45, 90]):
            lampam_weightings[3] = 0
            lampam_weightings[7] = 0
            lampam_weightings[11] = 0

        # modifying lamination parameter factor for orthotropy requirements
        if constraints.bal:
            lampam_weightings[2] = 0
            lampam_weightings[3] = 0
        if constraints.oopo:
            lampam_weightings[10] = 0
            lampam_weightings[11] = 0

        mean = np.average(lampam_weightings[lampam_weightings != 0])

        if constraints.oopo:
            if set(constraints.set_of_angles) == set([0, 45, -45, 90]):
                lampam_weightings[10] = mean * obj_func_param.coeff_oopo
            else:
                lampam_weightings[10] = mean * obj_func_param.coeff_oopo
                lampam_weightings[11] = mean * obj_func_param.coeff_oopo

#        if constraints.bal and obj_func_param.penalty_ipo_switch
#        and (obj_func_param.penalty_bal_ipo_switch_mp or
#             (not obj_func_param.penalty_bal_ipo_switch_mp and is_thick)):
#
#        lampam_weightings_3 = np.copy(lampam_weightings)
#
#        if constraints.bal and obj_func_param.penalty_ipo_switch:
#            if set(constraints.set_of_angles) == set([0, 45, -45, 90]):
#                lampam_weightings_3[2] = mean * obj_func_param.coeff_bal_ipo
#            else:
#                lampam_weightings_3[2] = mean * obj_func_param.coeff_bal_ipo
#                lampam_weightings_3[3] = mean * obj_func_param.coeff_bal_ipo
#
#        if constraints.oopo:
#            if set(constraints.set_of_angles) == set([0, 45, -45, 90]):
#                lampam_weightings_3[10] = mean * obj_func_param.coeff_oopo
#            else:
#                lampam_weightings_3[10] = mean * obj_func_param.coeff_oopo
#                lampam_weightings_3[11] = mean * obj_func_param.coeff_oopo
#
#        self.sum_lampam_weightings_3 = np.sum(lampam_weightings_3)
#        self.lampam_weightings_3 \
#        = lampam_weightings_3 / self.sum_lampam_weightings_3

#        if not np.allclose(lampam_weightings, self.lampam_weightings_ini):
#            print(f"""
#The lamination-parameter weightings have been modified (before normalisation):
#{self.lampam_weightings_ini} -> {lampam_weightings}
#    """)

        self.sum_lampam_weightings = np.sum(lampam_weightings)
        self.lampam_weightings = lampam_weightings \
        / self.sum_lampam_weightings

        self.lampam_weightingsA = self.lampam_weightings[0:4]
        self.lampam_weightingsD = self.lampam_weightings[8:12]

    def set_n_plies(self, n_plies, constraints):
        """
        returns the number of plies of a laminate
        and the number of a its potential middle ply (0 otherwise)
        """

        if not isinstance(n_plies, int):
            raise PanelDefinitionError("""
Attention, the number of plies in the panel must be an integer!""")
        if n_plies < 1:
            raise PanelDefinitionError("""
Attention, the number of plies in the panel must be positive!""")
        middle_ply_index = 0
        self.has_middle_ply = False
        if constraints.sym:
            if n_plies % 2 == 1:
                middle_ply_index = int((n_plies + 1)/2)
                self.has_middle_ply = True
        self.n_plies, self.middle_ply_index = n_plies, middle_ply_index

#        if constraints.sym and self.middle_ply_index != 0:
#            raise PanelDefinitionError("""
#Attention, the number of plies in the panel must be even for
#guide-based-blending !""")

        return 0

    def calc_weight(self, density_area):
        """
returns the weight of a panel
        """
        return density_area*self.area*self.n_plies


    def show(self):
        " Display object - non verbose verisn or __repr__"
        return f"""
Panel ID : {self.ID}
Number of plies : {self.n_plies}
"""


    def __repr__(self):
        " Display object "
        to_disp = ''
        # Lamination-parameter targets
        if np.array(self.lampam_target).size:
            to_disp = to_disp \
            + 'Lamination-parameter targets : ' \
            + str(self.lampam_target)  + '\n'

        # Lamination-parameter weighting in the panel objective function
        if hasattr(self, 'lampam_weightings2'):
            to_disp = to_disp \
            + """
Final lamination-parameter weighting in the panel objective function: :
    A : {self.lampam_weightings2[0:4]}
    B : {self.lampam_weightings2[4:8]}
    D : {self.lampam_weightings2[8:12]}
    """

        return f"""
Panel ID : {self.ID}
Number of plies : {self.n_plies}
Neighbour panel IDs: {self.neighbour_panels}
Length in the x-direction : {self.length_x}
Length in the y-direction : {self.length_y}
Area : {self.area}
Weighting in multi-panel objective funcion : {self.weighting}
Load intensity in th x-direction : {self.N_x}
Load intensity in th y-direction : {self.N_y}
Position of potential middle ply : {self.middle_ply_index}
    """ + to_disp

#Lamination-parameter weighting in the panel objective function:
#    during blending steps 3 and 4.1
#        A : {self.lampam_weightings[0:4]}
#        B : {self.lampam_weightings[4:8]}
#        D : {self.lampam_weightings[8:12]}
#    during blending steps 4.2 and 4.3
#        A : {self.lampam_weightings2[0:4]}
#        B : {self.lampam_weightings2[4:8]}
#        D : {self.lampam_weightings2[8:12]}


class PanelDefinitionError(Exception):
    " Errors during the definition of a panel"

if __name__ == "__main__":
    print('*** Test for the class Panel ***\n')
    constraints = Constraints(sym=True)
    panel1 = Panel(1, constraints, n_plies=12)
    print(panel1)
