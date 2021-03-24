# -*- coding: utf-8 -*-
"""
Function to check the feasibility of laminate lay-ups for the contiguity design
guideline.
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

def internal_contig2(stack, constraints):
    '''
    returns True if a laminate lay-up staisfy the contiguity design guideline,
    False otherwise

    INPUTS

    - stack: the laminate stacking sequence
    - constraints: the set of constraints
    '''
    if stack.ndim == 2:
        stack = stack.reshape((stack.size, ))

    if not constraints.contig:
        return True

    if constraints.n_contig < stack.size:
        diff = stack.size - constraints.n_contig

        if constraints.n_contig == 2:
            for jj in np.arange(diff):
                if stack[jj]==stack[jj + 1] \
                and stack[jj]==stack[jj + 2]:
                    return False

        elif constraints.n_contig == 3:
            for jj in np.arange(diff):
                if stack[jj]==stack[jj + 1] \
                and stack[jj]==stack[jj + 2] \
                and stack[jj]==stack[jj + 3]:
                    return False

        elif constraints.n_contig == 4:
            for jj in np.arange(diff):
                if stack[jj]==stack[jj + 1] \
                and stack[jj]==stack[jj + 2] \
                and stack[jj]==stack[jj + 3] \
                and stack[jj]==stack[jj + 4]:
                    return False

        elif constraints.n_contig == 5:
            for jj in np.arange(diff):
                if stack[jj]==stack[jj + 1] \
                and stack[jj]==stack[jj + 2] \
                and stack[jj]==stack[jj + 3] \
                and stack[jj]==stack[jj + 4] \
                and stack[jj]==stack[jj + 5]:
                    return False

        elif constraints.n_contig == 6:
            for jj in np.arange(diff):
                if stack[jj]==stack[jj + 1] \
                and stack[jj]==stack[jj + 2] \
                and stack[jj]==stack[jj + 3] \
                and stack[jj]==stack[jj + 4] \
                and stack[jj]==stack[jj + 5] \
                and stack[jj]==stack[jj + 6]:
                    return False
        else:
            raise Exception(
                    'constraints.n_contig must be 2, 3, 4 or 5')
    return True


def internal_contig(angle, constraints, angle2=None):
    '''
returns only the stacking sequences that satisfy the contiguity rule

OUTPUTS

- angle: the selected sublaminate stacking sequences line by
line
- angle2: the selected sublaminate stacking sequences line by
line if a second sublaminate is given as input for angle2

INPUTS

- angle: the first sublaminate stacking sequences
- angle:2 matrix storing the second sublaminate stacking sequences

    '''
    if angle.ndim == 1:
        angle = angle.reshape((1, angle.size))
    n_plies_group = angle.shape[1]

    # TO ENSURE CONTIGUITY
    if constraints.contig:
        # To ensure the contiguity constraint within groups of plies

        if not angle2 is None:

            if constraints.n_contig < n_plies_group:
                diff = n_plies_group-constraints.n_contig

                if constraints.n_contig == 2:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj + 2]==angle[ii, jj + 1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break

                elif constraints.n_contig == 3:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj + 2]==angle[ii, jj + 1] \
                            and angle[ii, jj + 3]==angle[ii, jj + 1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break

                elif constraints.n_contig == 4:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj + 2]==angle[ii, jj + 1] \
                            and angle[ii, jj + 2]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break

                elif constraints.n_contig == 5:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4] \
                            and angle[ii, jj]==angle[ii, jj + 5]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break

                elif constraints.n_contig == 6:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4] \
                            and angle[ii, jj]==angle[ii, jj + 5] \
                            and angle[ii, jj]==angle[ii, jj + 6]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break

                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')


        else:

            if constraints.n_contig < n_plies_group:
                diff = n_plies_group-constraints.n_contig

                if constraints.n_contig == 2:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break

                elif constraints.n_contig == 3:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break

                elif constraints.n_contig == 4:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break

                elif constraints.n_contig == 5:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4] \
                            and angle[ii, jj]==angle[ii, jj + 5]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break

                elif constraints.n_contig == 6:

                    a = angle.shape[0]
                    for ii in range(a)[::-1]:

                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4] \
                            and angle[ii, jj]==angle[ii, jj + 5] \
                            and angle[ii, jj]==angle[ii, jj + 6]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break

                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')

    return angle, angle2


if __name__ == "__main__":
    'Test'

    import sys
    sys.path.append(r'C:\BELLA')
    from src.LAYLA_V02.constraints import Constraints
    from src.divers.pretty_print import print_list_ss

    constraints = Constraints()
    constraints.contig = True
    constraints.n_contig = 5

    print('*** Test for the function internal_contig ***\n')
    print('Input stacking sequences:\n')
    ss = np.array([[0, 0, 45, 0, 0, 45, 90, 0, 45, 90],
                   [0, 0, 0, 0, 0, 0, 90, 90, 45, 90],
                   [0, 0, 0, 0, 0, 0, 90, 90, 45, 90]])
    print_list_ss(ss)
    test = internal_contig(ss, constraints)[0]
    if test.shape[0]:
        print('Stacking sequences satisfying the rule:\n')
        print_list_ss(test)
    else:
        print('No stacking sequence satisfy the rule\n')




