# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 20:03:35 2020

@author: slaiad
"""
from ProcessDbPack.ParsePep import getAllProts, getAllPepsForParallel
from ProcessSpecPack.ParseMsms import readPeaksFromXML, readTagsFromMS2
from ProcessSpecPack.PreProcess import deIsotope
import concurrent.futures as cf
from itertools import repeat
import pandas as pd
import numpy as np
from Utils.Funcs import binarySearch
from ScorePack.Match import getPepCand
from Utils.Funcs import divideInputList

#%% Set Parameters
from Parameters import *
#%% functions

def getFeasiblePep(pcMass, allPeps, massTol):
    lb = binarySearch(allPeps['pcMass'], 0, len(allPeps)-1, pcMass*(1 - massTol))
    ub = binarySearch(allPeps['pcMass'], 0, len(allPeps)-1, pcMass*(1 + massTol))
    # print(lb, ub)
    return list(allPeps['seq'][lb : ub])
    
# def searchOne(scanNo):
def searchSome(scanNoList):
    global specDict, allPeps
    psms = []
    # i = 0
    for scanNo in scanNoList:
        # i+=1
        # # noneRes = tuple((scanNo, [], []))
        # # global specDict, allPeps
        # if i % 200 == 0:
        #     print('processing i files ', i, scanNo)
            
        pcMass = specDict[scanNo][0]
        rawSpec = specDict[scanNo][1]
        mzs = deIsotope(rawSpec)
        
        if len(mzs) == 0:
            continue
        
        # print(mzs)
        
        tags = readTagsFromMS2(mzs)
        
        # print('finish tag', scanNo)
            
        reliableTags = []
        for tag in tags:
            
            if tag[1] < TAG_SCORE_THRES:
                break
            
            reliableTags.append(tag[0])
                
            if len(reliableTags) == MAX_TAGS_NUM:
                break
            
        if len(reliableTags) == 0:
            continue
        # print('finish reliableTags', scanNo)
        
        reliableTags = list(set(reliableTags))
        
        feasiblePeps = getFeasiblePep(pcMass, allPeps, MASS_TOL)
        # print(scanNo, pcMass)
        # print('reliable tags')
        # print(reliableTags)
        # print('length', len(reliableTags), len(set(reliableTags)))
        # print('IAHYNKR in feasible', 'IAHYNKR' in feasiblePeps)
        pepCand = getPepCand(reliableTags, feasiblePeps)
        # print('pepCand', pepCand)
        # print('IAHYNKR in pepCand', 'IAHYNKR' in pepCand)
        psms.append(tuple((scanNo, reliableTags, pepCand, len(feasiblePeps))))  
        # print('MYSYPARVPPPPPIAR' in pepCand)
    return psms
    
    
def doSearch(specDict, allPeps):
    
    # tagDict = {}
    # pepCandDict = {}
    scanTagCand = {}
    
    scanNoList = [key for key in specDict]
    
    # scanNoList = [26477] #test
    
    dividedScanList = divideInputList(scanNoList, 2000)
    
    print('len of scanNo list ', len(scanNoList))
    allPsms = []
    print("lsz start computation")
    with cf.ProcessPoolExecutor(max_workers = 88) as executor:
        # results = executor.map(searchOne, scanNoList)
        for result in executor.map(searchSome, dividedScanList):
        # for result in results:
            # if len(result[2]) != 0:
            allPsms.extend(result)
                # print(result)
                # break
                # scanTagCand[result[0]] = [result[1], result[2]]
    print('lsz finished computation')

    return allPsms

def parseDB(protList):
    dividedList = divideInputList(protList, 100)
    
    dtype = [('seq', np.unicode_, MAX_LENGTH_OF_PEP), ('pcMass', float)]
    allPeps = np.array([], dtype)
    
    with cf.ProcessPoolExecutor(max_workers = 88) as executor:
        for result in executor.map(getAllPepsForParallel, dividedList):
            allPeps = np.concatenate((allPeps, result), axis = 0)

    allPeps = np.unique(allPeps)        
    allPeps = (np.sort(allPeps, order = ['pcMass']))
    print("Number of all peps = ", len(allPeps))
    
    return allPeps

#%% Read DB and Spectra
if __name__ == '__main__':
    from time import time
    
    protList = getAllProts(DB_PATH)
    # allPeps = getAllPeps(protList)
    t2 = time()
    
    
    allPeps = parseDB(protList)
    
    t3 = time()        
    print("finishe parsing database ", t3-t2)
    print("Number of all peps = ", len(allPeps))
    
    specDict = readPeaksFromXML(SPEC_PATH)
    
    #%% doSearch
    from Parameters import *
    print("start Searching")
    
    t4 = time()
    allPsms = doSearch(specDict, allPeps)
    t5 = time()
    
    print("finish searching ", t5-t4)
    from Utils.Funcs import sendEmail
    # sendEmail()
    #% write
    dataPsms = pd.DataFrame({'scanNo' : [psm[0] for psm in allPsms],
                                    'tags'   : [psm[1] for psm in allPsms],
                                    'seqs'   : [psm[2] for psm in allPsms]})#, orient='scanNo',columns=['Tags', 'Candidates'])
    filePath = pd.ExcelWriter(r'TempData/dataPsms.xlsx')
    dataPsms.to_excel(r'TempData/dataPsms.xlsx')
    
    #% verification
    dataMascot = pd.read_excel('/home/slaiad/Code/TagTree/testData/MK_SIO13_P2_GM1.xlsx')
    
    masDict = {}
    # cometModDict = {}
    for i in range(len(dataMascot['scanNo'])):
        pep = dataMascot['pep'][i]
        masDict[dataMascot['scanNo'][i]] = pep
        
    result = []
    numSameMascot = 0
    numDiffMascot = 0
    numNoInMasCot = 0
    for psm in allPsms:
        pepCand = psm[2]
        scanNo = psm[0]
        if scanNo in masDict:
            pep = masDict[scanNo].replace("L", "I")
            if pep.upper() in [cand.upper() for cand in pepCand]:
                numSameMascot += 1
            else:
                numDiffMascot += 1
                # print(scanNo, masDict[scanNo], pepCand)
            result.append([scanNo, pepCand, psm[1], masDict[scanNo], pep.upper() in [cand.upper() for cand in pepCand], psm[3]])
        else:
            numNoInMasCot += 1
    print('total number of res same as Mascot', numSameMascot)
    print('total number of res diff from Mascot', numDiffMascot)
    print('total number of res not in Mascot', numNoInMasCot)
    
    #%
    dataVeri = pd.DataFrame({'scanNo'     : [res[0] for res in result],
                             'pepCand'    : [res[1] for res in result],
                             'tags'       : [res[2] for res in result],
                             'truth'      : [res[3] for res in result],
                             'inOrnot'    : [res[4] for res in result],
                             'numFesiPep' : [res[5] for res in result]
                             })#, orient='scanNo',columns=['Tags', 'Candidates'])
    # filePath = pd.ExcelWriter(r'TempData/dataVeri.xlsx')
    dataVeri.to_excel(r'TempData/tag1p6_AA0p5_pcTOL0.1_mixFuzz80top20set_QeGA_bigDB_nTerm42p011_varOxiM15p99.xlsx')
    
    
    
   

