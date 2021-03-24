# -*- coding: utf-8 -*-
"""
Class for the parameters of the optimiser
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np


class Parameters():
    " An object for storing the optimiser parameters "

    def __init__(
            self,
            constraints,
            group_size_min=4,
            group_size_max=12,
            global_node_limit=100,
            local_node_limit=100,
            local_node_limit_final=1,
            global_node_limit_final=5,
            global_node_limit2=10,
            local_node_limit2=10,
            global_node_limit3=10,
            local_node_limit3=10,
            n_ini_ply_drops=5,
            p_A=100,
            n_D1=6,
            n_D2=10,
            n_D3=1,
            repair_membrane_switch=True,
            repair_flexural_switch=True,
            n_plies_ref_panel=1000,
            save_success_rate=False,
            time_limit_group_pdl=1,
            time_limit_all_pdls=100,
            save_buckling=False):
        " Create a set of parameters for BELLA"

        self.save_buckling = save_buckling

        ### BELLA step 2

        # size of the ply groups during BELLA step 2
        #   desired group size for smaller groups
        self.group_size_min = np.around(group_size_min)
        if not isinstance(group_size_min, int):
            raise ParametersDefinitionError("""
Attention, group_size_min must be an integer!""")
        if group_size_min < 1:
            raise ParametersDefinitionError("""
Attention, group_size_min must be strictly positive!""")
        #   maximum group size
        if not isinstance(group_size_max, int):
            raise ParametersDefinitionError("""
Attention, group_size_max must be an integer!""")
        self.group_size_max = group_size_max
        if group_size_min > group_size_max:
            raise ParametersDefinitionError("""
Attention, group_size_min must smaller than group_size_max!""")

        # Time limit to create a group ply-drop layout
        self.time_limit_group_pdl = time_limit_group_pdl
        # Time limit to create a ply-drop layout
        self.time_limit_all_pdls = time_limit_all_pdls

        # Number of initial ply drops to be tested
        self.n_ini_ply_drops = np.around(n_ini_ply_drops)
        if not isinstance(n_ini_ply_drops, int):
            raise ParametersDefinitionError("""
Attention, n_ini_ply_drops must be an integer!""")
        if n_ini_ply_drops < 1:
            raise ParametersDefinitionError("""
Attention, n_ini_ply_drops must be strictly positive!""")


        ### BELLA step 3

        # Branching limits for global pruning during BELLA step 3
        self.global_node_limit = np.around(global_node_limit)
        if not isinstance(global_node_limit, int):
            raise ParametersDefinitionError("""
Attention, global_node_limit must be an integer!""")
        if global_node_limit < 1:
            raise ParametersDefinitionError("""
Attention, global_node_limit must be strictly positive!""")
        self.global_node_limit_final = np.around(global_node_limit_final)
        if not isinstance(global_node_limit_final, int):
            raise ParametersDefinitionError("""
Attention, global_node_limit_final must be an integer!""")
        if global_node_limit_final < 1:
            raise ParametersDefinitionError("""
Attention, global_node_limit_final must be strictly positive!""")
        self.local_node_limit = np.around(local_node_limit)
        if not isinstance(local_node_limit, int):
            raise ParametersDefinitionError("""
Attention, local_node_limit must be an integer!""")
        if local_node_limit < 1:
            raise ParametersDefinitionError("""
Attention, local_node_limit must be strictly positive!""")
        self.local_node_limit_final = np.around(
        local_node_limit_final)
        if not isinstance(local_node_limit_final, int):
            raise ParametersDefinitionError("""
Attention, local_node_limit_final must be an integer!""")
        if local_node_limit_final < 1:
            raise ParametersDefinitionError("""
Attention, local_node_limit_final must be strictly positive!""")







        ### BELLA step 4.1

        ### Thickness of the reference panels
        if isinstance(n_plies_ref_panel, (int, float)):
            self.n_plies_ref_panel = n_plies_ref_panel
        else:
            raise ParametersDefinitionError("""
The ply count of the reference panels must be a number!""")

        # repair to improve the convergence of in-plane lamination parameters
        # and of out-of-plane lamination parameters
        self.repair_membrane_switch = repair_membrane_switch
        if not isinstance(repair_membrane_switch, bool):
            raise ParametersDefinitionError("""
Attention, repair_membrane_switch should be a boolean value!""")
        self.repair_flexural_switch = repair_flexural_switch
        if not isinstance(repair_flexural_switch, bool):
            raise ParametersDefinitionError("""
Attention, repair_flexural_switch should be a boolean value!""")

        self.save_success_rate = save_success_rate

        # coefficient for the proportion of the laminate thickness that can be
        # modified during the refinement for membrane properties in the repair
        # process
        self.p_A = p_A
        if not isinstance(p_A, (float, int)):
            raise ParametersDefinitionError("""
Attention, p_A should have a numeric value!""")
        if not (0 <= p_A <= 100):
            raise ParametersDefinitionError("""
Attention, p_A must be between 0 and 100!""")

        # n_D1: number of plies in the last permutation
        # during repair for disorientation and/or contiguity
        self.n_D1 = n_D1
        if not isinstance(n_D1, int):
            raise ParametersDefinitionError("""
Attention, n_D1 must be an integer!""")

        # n_D2: number of ply shifts tested at each step of the
        # re-designing process during refinement of flexural properties
        self.n_D2 = n_D2
        if not isinstance(n_D2, int):
            raise ParametersDefinitionError("""
Attention, n_D2 must be an integer!""")

        # n_D3: number of times the algorithms 1 and 2 are repeated during the
        # flexural property refinement
        self.n_D3 = n_D3
        if not isinstance(n_D3, int):
            raise ParametersDefinitionError("""
Attention, n_D2 must be an integer!""")


        ### BELLA step 4.2

        # Branching limits for global pruning during BELLA step 4.2
        self.global_node_limit2 = np.around(global_node_limit2)
        if not isinstance(global_node_limit2, int):
            raise ParametersDefinitionError("""
Attention, global_node_limit2 must be an integer!""")
        if global_node_limit2 < 1:
            raise ParametersDefinitionError("""
Attention, global_node_limit2 must be strictly positive!""")
        self.local_node_limit2 = np.around(local_node_limit2)
        if not isinstance(local_node_limit2, int):
            raise ParametersDefinitionError("""
Attention, local_node_limit2 must be an integer!""")
        if local_node_limit2 < 1:
            raise ParametersDefinitionError("""
Attention, local_node_limit2 must be strictly positive!""")


        ### BELLA step 4.3

        # Branching limits for global pruning during BELLA step 4.3
        self.global_node_limit3 = np.around(global_node_limit3)
        if not isinstance(global_node_limit3, int):
            raise ParametersDefinitionError("""
Attention, global_node_limit3 must be an integer!""")
        if global_node_limit3 < 1:
            raise ParametersDefinitionError("""
Attention, global_node_limit3 must be strictly positive!""")
        self.local_node_limit3 = np.around(local_node_limit3)
        if not isinstance(local_node_limit3, int):
            raise ParametersDefinitionError("""
Attention, local_node_limit3 must be an integer!""")
        if local_node_limit3 < 1:
            raise ParametersDefinitionError("""
Attention, local_node_limit3 must be strictly positive!""")


    def __repr__(self):
        " Display object "

        return f"""
Parameters of BELLA step 2:
    Number of initial ply drops to be tested: {self.n_ini_ply_drops}
    Minimum size for the ply groups: {self.group_size_min}
    Maximum size for the ply groups: {self.group_size_max}
    Time limit to generate pdl at each group search: {self.time_limit_group_pdl}
    Time limit to generateall initial pdls: {self.time_limit_all_pdls}

Parameters of BELLA step 3:
    Branching limits during beam search:
        - for global pruning: {self.global_node_limit}
        - for local pruning: {self.local_node_limit}
        - for global pruning at the last level: {self.global_node_limit_final}
        - for local pruning at the last level: {self.local_node_limit_final}

Parameters of BELLA step 4.1:
    Input number of plies in reference panel: {self.n_plies_ref_panel}
    Repair for in-plane lamination parameter convergence: {self.repair_membrane_switch}
    Repair for out-of-plane lamination parameter convergence: {self.repair_flexural_switch}
    p_A: {self.p_A}%
    n_D1: {self.n_D1}
    n_D2: {self.n_D2}
    n_D3: {self.n_D3}

Parameters of BELLA step 4.2:
    Branching limits during beam search:
        - for global pruning: {self.global_node_limit2}
        - for local pruning: {self.local_node_limit2}

Parameters of BELLA step 4.3:
    Branching limits during beam search:
        - for global pruning: {self.global_node_limit3}
        - for local pruning: {self.local_node_limit3}
"""

class ParametersDefinitionError(Exception):
    """ Error during parameter definition"""

if __name__ == "__main__":
    import sys
    sys.path.append(r'C:\BELLA')
    from src.BELLA.constraints import Constraints
    constraints = Constraints(sym=True)
    constraints.bal = True
    constraints.oopo = True
    parameters = Parameters(
        constraints=constraints)
    print(parameters)
