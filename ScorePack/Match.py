#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 09:53:32 2021

@author: slaiad@ust.hk
"""

from fuzzywuzzy import fuzz
import numpy as np
from Parameters import MAX_LENGTH_OF_PEP

def simi(str1, str2):
    # print(max(fuzz.partial_ratio(str1, str2) , 
    #            fuzz.partial_ratio(str1[::-1], str2)),'+',len(str1)/len(str2)*10)
    # return max(fuzz.partial_ratio(str1, str2) + fuzz.ratio(str1, str2)/8, 
    #         fuzz.partial_ratio(str1[::-1], str2) + fuzz.ratio(str1[::-1], str2)/8)
    
    if str1 in str2 or str1[::-1] in str2:
        return 100 - min((len(str2) - len(str1)), 15)
    return max(fuzz.partial_ratio(str1, str2) , 
               fuzz.partial_ratio(str1[::-1], str2))+ len(str1)/len(str2)*10


def getPepCand(reliableTags, feasiblePeps):
    
    dtype = [('seq', np.unicode_, MAX_LENGTH_OF_PEP), ('simi', float)]
    pepCands = []
    
    additionalTags = []
    
    for tag in reliableTags:
        if tag.count('Q') == 1:
            additionalTags.append(tag.replace('Q','GA'))
            additionalTags.append(tag.replace('Q','AG'))
    
    reliableTags.extend(additionalTags)
    for tag in reliableTags:
            
        # if tag.count('GA')
        
        simiFuzz = np.array([], dtype)
        
        for pep in feasiblePeps:
            # if tag[0] in pep or tag[0][::-1] in pep:# exact containing
            similarity = simi(tag, pep)
            # if pep == 'MYSYPARVPPPPPIAR':
            #     print(tag, similarity)
            # if pep == 'IAHYNKR':
            #     print(tag, similarity)
            if similarity >= 80:
                simiFuzz = np.concatenate((simiFuzz, np.array([(pep, similarity)], dtype)), axis = 0)
        sortSimiFuzz = np.sort(simiFuzz, order = 'simi')
        # print(sortSimiFuzz)
        pepCands.extend(list(sortSimiFuzz['seq'][ : - min(20, len(sortSimiFuzz) + 1) : -1]))
        # print(simiFuzz)       
    # print(pepCands)
    return list(set(pepCands))
    # return list(sortSimiFuzz['seq'][ - min(20, len(sortSimiFuzz)) : ])


if __name__ == '__main__':
    
    # print(simi('QIWLEJSBS', 'IQWLEJSBS'))
    
    s1 = 'AAE'
    
    s2 = 'ALAEEAAKKGR'
    print(simi(s1,s2))
    
    print(fuzz.partial_ratio(s1, s2))