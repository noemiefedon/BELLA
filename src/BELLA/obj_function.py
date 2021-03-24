# -*- coding: utf-8 -*-
"""
Class for the objective function parameters
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np


class ObjFunction():
    " An object for storing the objective function parameters"

    def __init__(
            self,
            constraints,
            coeff_diso=0,
            coeff_contig=0,
            coeff_10=0,
            coeff_bal_ipo=1,
            coeff_oopo=1,
            coeff_spacing=1,
            n_modes=1):
        "Initialise the objective function parameters"

        # number of buckling modes to be tested
        self.n_modes = n_modes
        if not isinstance(n_modes, int):
            raise ObjFunctionDefinitionError("""
Attention, n_modes must be an integer!""")
        if n_modes < 1:
            raise ObjFunctionDefinitionError("""
Attention, n_modes must be strictly positive!""")

        ### 10% rule
        if not constraints.rule_10_percent:
            self.coeff_10 = 0
        else:
            self.coeff_10 = coeff_10
            if not isinstance(coeff_10, (float, int)):
                raise ObjFunctionDefinitionError("""
Attention, coeff_10 must be a number (float or integer)!""")
            if coeff_10 < 0:
                raise ObjFunctionDefinitionError("""
The weight of penalty for the 10% rule must be a positive!""")

        ### balance
        if not constraints.bal:
            self.coeff_bal_ipo = 0
        else:
            self.coeff_bal_ipo = coeff_bal_ipo
            if not isinstance(coeff_bal_ipo, (float, int)):
                raise ObjFunctionDefinitionError("""
Attention, coeff_bal_ipo must be a number (float or integer)!""")
            if coeff_bal_ipo < 0:
                raise ObjFunctionDefinitionError("""
The weight of penalty for in-plane orthotropy must be a positive!""")

        ### out-of-plane orthotropy
        if not constraints.oopo:
            self.coeff_oopo = 0
        else:
            self.coeff_oopo = coeff_oopo
            if not isinstance(coeff_oopo, (float, int)):
                raise ObjFunctionDefinitionError("""
Attention, coeff_oopo must be a number (float or integer)!""")
            if coeff_oopo < 0:
                raise ObjFunctionDefinitionError("""
The weight of penalty for out-of-plane orthotropy must be a positive!""")

        ### contiguity rule
        if not constraints.contig:
            self.coeff_contig = 0
        else:
            self.coeff_contig = coeff_contig
            if coeff_contig < 0:
                raise ObjFunctionDefinitionError("""
The weight of the contiguity constraint penalty must be a positive!""")

        ### disorientation rule
        if not constraints.diso:
            self.coeff_diso = 0
        else:
            self.coeff_diso = coeff_diso
            if coeff_diso < 0:
                raise ObjFunctionDefinitionError("""
The weight of the disorientation constraint penalty must be a positive!""")


        ### ply drop spacing guideline
        if not constraints.pdl_spacing:
            self.coeff_spacing = 0
        else:
            self.coeff_spacing = coeff_spacing
            if not isinstance(coeff_spacing, (float, int)):
                raise ObjFunctionDefinitionError("""
The weight of penalty for the ply drop spacing guideline must be a number!""")
            if coeff_spacing < 0:
                raise ObjFunctionDefinitionError("""
The weight of penalty for the ply drop spacing guideline must be positive!""")

#        ### ply-drop guidelines
#        ### weight of the panel boundaries when calculating ply-drop layout
#        # penalties
#        # 1 for boundaries of level 1
#        # 1 - coeff for boundaries of level 2
#        # 1 - 2 * coeff for boundaries of level 3 ...
#        if not constraints.pdl_spacing:
#            self.coeff_panel_pdls = 0
#        else:
#            self.coeff_panel_pdls = coeff_panel_pdls


    def set_initial_panel_weightings(self, multipanel):
        """
        to set the initial panel weightings
        """
        self.reduced_panel_weightings \
        = np.array([p.weighting for p in multipanel.reduced.panels])

        self.panel_weightings_ini \
        = np.array([p.weighting for p in multipanel.panels])

    def __repr__(self):
        " Display object "

        return f"""
Number of buckling modes to be tested if buckling minimisation problem: {self.n_modes}

Penalty coefficients:
    - for the disorientation rule: {self.coeff_diso}
    - for the contiguity rule: {self.coeff_contig}
    - for the 10% rule: {self.coeff_10}
    - for balance: {self.coeff_bal_ipo}
    - for out-of-plane orthotropy: {self.coeff_oopo}
    - for the ply drop spacing guideline: {self.coeff_spacing}
"""
#    - for the weight of the panel boundaries when calculating ply-drop layout
#    penalties: {self.coeff_panel_pdls}


class ObjFunctionDefinitionError(Exception):
    """ Error during parameter definition"""

if __name__ == "__main__":
    import sys
    sys.path.append(r'C:\BELLA')
    from src.BELLA.constraints import Constraints
    constraints = Constraints(sym=True)
    constraints.bal = True
    constraints.oopo = True
    obj_func_param = ObjFunction(
        constraints=constraints,
        lampam_weightings=np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]),
        coeff_bal_ipo=1000,
        coeff_oopo=20)
    print(obj_func_param)
