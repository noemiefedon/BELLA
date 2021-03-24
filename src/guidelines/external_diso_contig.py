# -*- coding: utf-8 -*-
"""
Function to ensure external sorting for disorientation and contiguity
'across adjacent sublaminates'

Created on Mon Jan 29 12:00:18 2018

@author: Noemie Fedon
"""
import numpy as np

def external_diso_contig(angle, n_plies_group, constraints, ss_before, angle2 = None):
    '''
returns only the stacking sequences that satisfy constraints concerning
contiguity and disorientation at the junction with an adjacent group
of plies, but not within the group of plies 

OUTPUTS

- angle: the selected sublaminate stacking sequences line by 
line
- angle2: the selected sublaminate stacking sequences line by 
line if a second sublaminate is given as input for angle2
 
INPUTS

- angle: the first sublaminate stacking sequences 
- angle:2 matrix storing the second sublaminate stacking sequences 
- n_plies_group: number of plies in the sublaminate (scalar)
- ss_before is the stacking sequence of the sublaminate adjacent to the first 
sublaminate

    '''
    if angle.ndim == 1:
        angle = angle.reshape((1, angle.size))
        
    ss_beforeLength = ss_before.size
    delta = constraints.delta_angle
    
    # CHECK FOR CORRECT INPUTS SIZE
    if n_plies_group > angle.shape[1]:
        raise Exception('The input set of angles have fewer elements that what is asked to be checked')

    if angle2 is None:
        # TO ENSURE DISORIENTATION
        if constraints.diso:
            # To ensure the disorientation constraint at the 
            # junction of groups of plies with precedent plies
            if ss_beforeLength != 0:
                a = angle.shape[0]
                for ii in range(a)[::-1]:
                    if abs(angle[ii, 0]-ss_before[-1])> delta \
                    and abs(180 +angle[ii, 0]-ss_before[-1])>delta  \
                    and abs(angle[ii, 0]-ss_before[-1]-180)>delta:
                        angle = np.delete(angle, np.s_[ii], axis=0)
                        continue
        
        # TO ENSURE CONTIGUITY
        if constraints.contig:
            
            # To ensure the contiguity constraint at the junction of ply groups 
            if ss_beforeLength>=1:
                
                if constraints.n_contig ==2:
           
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if  angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
                            
                elif constraints.n_contig == 3:
                    
                    if n_plies_group>2:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
    
                elif constraints.n_contig == 4:
                    
                    if n_plies_group>3:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue      
                    
                elif constraints.n_contig == 5:
                    
                    if n_plies_group>4:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1]\
                            and angle[ii, 4] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
                            
                elif constraints.n_contig == 6:
                    
                    if n_plies_group>5:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1] \
                            and angle[ii, 4] == ss_before[-1] \
                            and angle[ii, 5] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
                    
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
            
    
            if ss_beforeLength>=2:
                if constraints.n_contig ==2:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if  angle[ii, 0] == ss_before[-1] \
                        and angle[ii, 0] == ss_before[-2]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            continue
                  
                elif constraints.n_contig == 3:
                    
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
               
                elif constraints.n_contig == 4:
                    
                    if n_plies_group>2:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
                  
                    
                elif constraints.n_contig == 5:
                    
                    if n_plies_group>3:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
                            
                elif constraints.n_contig == 6:
                    
                    if n_plies_group>4:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1] \
                            and angle[ii, 4] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
                            
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
            
    
            if ss_beforeLength>=3:
                if constraints.n_contig == 2:
                    pass
                elif constraints.n_contig == 3:
                    
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if angle[ii, 0] == ss_before[-1] \
                        and ss_before[-3] == ss_before[-1] \
                        and ss_before[-2] == ss_before[-1]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            continue
    
                elif constraints.n_contig == 4:
                    
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
    
                elif constraints.n_contig == 5:
                    
                    if n_plies_group>2:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
                            
                elif constraints.n_contig == 6:
                    
                    if n_plies_group>3:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1]  \
                            and angle[ii, 3] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
    
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
    
            if ss_beforeLength>=4:
                if constraints.n_contig ==2:
                    pass
                
                elif constraints.n_contig == 3:
                    pass
                
                elif constraints.n_contig == 4:
        
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if angle[ii, 0] == ss_before[-1] \
                        and ss_before[-4] == ss_before[-1] \
                        and ss_before[-3] == ss_before[-1] \
                        and ss_before[-2] == ss_before[-1]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            continue
    
                elif constraints.n_contig == 5:
                    
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and ss_before[-4] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
                            
                elif constraints.n_contig == 6:
                    
                    if n_plies_group>2:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and ss_before[-4] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue            
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
    
            if ss_beforeLength>=5:
                if constraints.n_contig ==2:
                    pass
                
                elif constraints.n_contig == 3:
                    pass
                
                elif constraints.n_contig == 4:
                    pass
                
                elif constraints.n_contig == 5:
                    
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if angle[ii, 0] == ss_before[-1] \
                        and ss_before[-4] == ss_before[-1] \
                        and ss_before[-3] == ss_before[-1] \
                        and ss_before[-2] == ss_before[-1] \
                        and angle[ii, 0] == ss_before[-5]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            continue
                            
                elif constraints.n_contig == 6:
                    
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and ss_before[-4] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 0] == ss_before[-5] \
                            and angle[ii, 1] == ss_before[-5]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                continue
    
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
            
            if ss_beforeLength>=6:
                if constraints.n_contig == 2:
                    pass
                
                elif constraints.n_contig == 3:
                    pass
                
                elif constraints.n_contig == 4:
                    pass
                
                elif constraints.n_contig == 5:
                    pass
 
                elif constraints.n_contig == 6:
                    
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if angle[ii, 0] == ss_before[-1] \
                        and ss_before[-4] == ss_before[-1] \
                        and ss_before[-3] == ss_before[-1] \
                        and ss_before[-2] == ss_before[-1] \
                        and angle[ii, 0] == ss_before[-5]  \
                        and ss_before[-6] == ss_before[-5]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            continue
                            
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')                    
        
    else:
    
        # TO ENSURE DISORIENTATION
        if constraints.diso:
            
            # To ensure the disorientation constraint at the 
            # junction of groups of plies with precedent plies
            if ss_beforeLength != 0:
                a = angle.shape[0]
                for ii in range(a)[::-1]:
                     
                    if abs(angle[ii, 0]-ss_before[-1])>delta \
                    and abs(180+angle[ii, 0]-ss_before[-1])>delta \
                    and abs(angle[ii, 0]-ss_before[-1]-180)>delta:
                        angle = np.delete(angle, np.s_[ii], axis=0)
                        angle2 = np.delete(angle2, np.s_[ii], axis=0)
                        continue
    
        # TO ENSURE CONTIGUITY
        if constraints.contig:
            
            # To ensure the contiguity constraint at the junction of ply groups
            if ss_beforeLength>=1:
                if constraints.n_contig == 2:
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if  angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
       
                elif constraints.n_contig == 3:
                    
                    if n_plies_group>2:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
      
                elif constraints.n_contig == 4:
                    
                    if n_plies_group>3:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
     
                elif constraints.n_contig == 5:
                    
                    if n_plies_group>4:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1] \
                            and angle[ii, 4] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
                
                elif constraints.n_contig == 6:
                    
                    if n_plies_group>5:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1] \
                            and angle[ii, 4] == ss_before[-1] \
                            and angle[ii, 5] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
              
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
    
            if ss_beforeLength>=2:
                if constraints.n_contig ==2:
                
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if  angle[ii, 0] == ss_before[-1] \
                        and angle[ii, 0] == ss_before[-2]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            angle2 = np.delete(angle2, np.s_[ii], axis=0)
                            continue
    
                elif constraints.n_contig == 3:
                    
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
    
                elif constraints.n_contig == 4:
                    
                    if n_plies_group>2:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
    
                elif constraints.n_contig == 5:
                    
                    if n_plies_group>3:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
     
                elif constraints.n_contig == 6:
                    
                    if n_plies_group>4:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1] \
                            and angle[ii, 4] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
        
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
    
            if ss_beforeLength>=3:
                if constraints.n_contig == 2:
                    pass
                        
                elif constraints.n_contig == 3:
                    
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if angle[ii, 0] == ss_before[-1] \
                        and ss_before[-3] == ss_before[-1] \
                        and ss_before[-2] == ss_before[-1]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            angle2 = np.delete(angle2, np.s_[ii], axis=0)
                            continue
    
                elif constraints.n_contig == 4:
                    
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
    
                elif constraints.n_contig == 5:
                    
                    if n_plies_group>2:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue

                elif constraints.n_contig == 6:
                    
                    if n_plies_group>3:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1] \
                            and angle[ii, 3] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
                            
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
    
            if ss_beforeLength >= 4:
                if constraints.n_contig == 2:
                    pass
                
                elif constraints.n_contig == 3:
                    pass
                
                elif constraints.n_contig == 4:
                    
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if angle[ii, 0] == ss_before[-1] \
                        and ss_before[-4] == ss_before[-1] \
                        and ss_before[-3] == ss_before[-1] \
                        and ss_before[-2] == ss_before[-1]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            angle2 = np.delete(angle2, np.s_[ii], axis=0)
                            continue
    
                elif constraints.n_contig == 5:
                    
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and ss_before[-4] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue

                elif constraints.n_contig == 6:
                    
                    if n_plies_group>2:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and ss_before[-4] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 1] == ss_before[-1] \
                            and angle[ii, 2] == ss_before[-1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
                            
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
    
            if ss_beforeLength>=5:
                if constraints.n_contig == 2:
                    pass
                
                elif constraints.n_contig == 3:
                    pass
                
                elif constraints.n_contig == 4:
                    pass
                
                elif constraints.n_contig == 5:
                    
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if angle[ii, 0] == ss_before[-1] \
                        and ss_before[-4] == ss_before[-1] \
                        and ss_before[-3] == ss_before[-1] \
                        and ss_before[-2] == ss_before[-1] \
                        and angle[ii, 0] == ss_before[-5]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            angle2 = np.delete(angle2, np.s_[ii], axis=0)
                            continue
 
                elif constraints.n_contig == 6:
                    
                    if n_plies_group>1:
                        a = angle.shape[0]
                        for ii in range(a)[::-1]:
                             
                            if angle[ii, 0] == ss_before[-1] \
                            and ss_before[-4] == ss_before[-1] \
                            and ss_before[-3] == ss_before[-1] \
                            and ss_before[-2] == ss_before[-1] \
                            and angle[ii, 0] == ss_before[-5]  \
                            and angle[ii, 1] == ss_before[-5]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                continue
                            
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')

            if ss_beforeLength>=6:
                if constraints.n_contig == 2:
                    pass
                
                elif constraints.n_contig == 3:
                    pass
                
                elif constraints.n_contig == 4:
                    pass
                
                elif constraints.n_contig == 5:
                    pass
 
                elif constraints.n_contig == 6:
                    
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                         
                        if angle[ii, 0] == ss_before[-1] \
                        and ss_before[-4] == ss_before[-1] \
                        and ss_before[-3] == ss_before[-1] \
                        and ss_before[-2] == ss_before[-1] \
                        and angle[ii, 0] == ss_before[-5]  \
                        and ss_before[-6] == ss_before[-5]:
                            angle = np.delete(angle, np.s_[ii], axis=0)
                            angle2 = np.delete(angle2, np.s_[ii], axis=0)
                            continue
                            
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
                    
    return angle, angle2
    
if __name__ == "__main__":  
    'Test'
    
    import sys
    sys.path.append(r'C:\BELLA')
    from src.LAYLA_V02.constraints import Constraints
    from src.divers.pretty_print import print_ss, print_list_ss
    
    constraints = Constraints()
    constraints.diso = True
    constraints.delta_angle = 45
    constraints.contig = True
    constraints.n_contig = 2

    print('*** Test for the function external_diso_contig ***\n')
    print('Input stacking sequences:\n')
    ss = np.array([[-45, -45, 0, 45, 90], 
                   [0, 45, 45, 45, 45],
                   [0, 0, 0, 45, 45]])
    print_list_ss(ss)
    print('Stacking sequence of adajacent sublaminate:\n')
    ss_before = np.array([-45])
    print_ss(ss_before)
    n_plies_group = 5 
    middle_ply = 0
    test, _ = external_diso_contig(ss, n_plies_group, constraints, ss_before, ss)
    if test.shape[0]:
        print('Stacking sequences satisfying the rule:\n')
        print_list_ss(test)
    else:
        print('No stacking sequence satisfy the rule\n')