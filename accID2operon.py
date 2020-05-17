import os
from pprint import pprint

import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO
import requests
import json
import time
import pickle

Entrez.email = os.getenv("EMAIL", "doelsnitz@utexas.edu")
headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 


def acc2MetaData(access_id: str):
    
    handle= Entrez.efetch(db='protein',id=access_id, rettype="ipg")
    
    #sometimes the 'ProteinList' Key is not present.
    #In that case, return a list with an empty dictionary
    try:
        proteinList = Entrez.read(handle)['IPGReport']["ProteinList"][0]
    except KeyError:
        print('ProteinList KeyError avoided')
        proteinList = [{}]

    protein = proteinList.get("CDSList", "MT")[0].attributes

    return [protein['accver'],protein['start'],protein['stop'],protein['strand']]

def getUID(NCacc, seq_start, seq_stop):
   
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term='+NCacc+'+AND+'+seq_start+':'+seq_stop+'[chrpos]&retmode=json', headers=headers)
    if response.ok: 
        uid = json.loads(response.text)["esearchresult"]["idlist"]
        try:
            uid = uid[0]
            uid = int(uid)
            return uid
        except:
            print("esearch couldn't find any associated uid")
    else:
        print('bad esearch request')


def uid2MetaData(uid):
    uid = str(uid)
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id='+uid+'&retmode=json', headers=headers)
    if response.ok: 
        data = json.loads(response.text)["result"][uid]
        description = data["description"]
        gdata = data["genomicinfo"][0]

        NCacc = gdata["chraccver"]
        seq_start = gdata["chrstart"]
        seq_stop = gdata["chrstop"]
        strand = getStrand(seq_start, seq_stop)
        return [description, NCacc, seq_start, seq_stop, strand, uid]
    else:
        print(response.status_code)


def getLinkID(uid):
    uid = str(uid)

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&id='+uid+'&retmode=json&idtype=acc', headers=headers)
    if response.ok: 
        P_ACC = json.loads(response.text)["linksets"][0]["linksetdbs"][0]["links"]
        try:
            P_ACC = P_ACC[0]
        except:
            print("elink couldn't find any accession associated with input uid")
            return
    else:
        print('bad elink request')
        return
    
    print(P_ACC)
    
    payload = {'from': 'P_REFSEQ_AC',
                'to': 'ACC',
                'format': 'list',
                'query': P_ACC,
                }

    response = requests.get('https://www.uniprot.org/uploadlists/', params=payload, headers=headers)
    if response.ok: 
        return response.text
    else:
        response.raise_for_status()



def getStrand(seq_start, seq_stop):
    if seq_stop > seq_start:
        return "+"
    elif seq_start > seq_stop:
        return "-"

def getTU(uid, seq_start, strand):

    '''
    Rules for inclusion/exclusion of genes from operon:
        - always take immediately adjacent genes
        - if query gene is in same direction as regulator, include it.
        - if query gene is expressed divergently from regulator, 
                grab all adjacent genes that are expressed divergently (change strand direction for next genes)
        - if query gene is expressed divergently from a co-transcribed neighbor of the regulaor, 
                grab that gene. (it may be another regulator. Important to know).
        - if query gene direction converges with regulator, exclude it.
    '''

    def getGene(geneStrand, direction, nextGene, geneList, UID):
        #print('starting '+direction)
        while geneStrand == nextGene[4]:
            #print(geneStrand, nextGene[4])
            if direction == '+':
                UID += 1
            elif direction == '-':
                UID -= 1
            nextGene = uid2MetaData(UID)
            if geneStrand == '-' and nextGene[4] == '+' and direction == '+':
                geneList.append(nextGene)
            elif geneStrand == nextGene[4]:
                geneList.append(nextGene)
            time.sleep(1)
    
    geneStrand = strand
    
    uidDOWN = uid-1
    downGene = uid2MetaData(uidDOWN)
    if seq_start > downGene[2]:
        #print('downGene is behind')
        if strand == '+' and downGene[4] == '-':
            geneStrand = downGene[4]
    else:
        pass
        #print('downGene is in front')
    downgenes = [downGene]
    getGene(geneStrand,'-',downGene, downgenes, uidDOWN)
    
    geneArray = list(reversed(downgenes))
    geneArray.append(uid2MetaData(uid))
    regulatorIndex = (len(geneArray)-1)

    geneStrand = strand

    uidUP = uid+1
    upGene = uid2MetaData(uidUP)
    if seq_start > upGene[2]:
        #print('upGene is behind')
        if strand == '+' and upGene[4] == '-':
            geneStrand = upGene[4]
    else:
        pass
        #print('upGene is in front')
    #upgenes = [upGene]
    geneArray.append(upGene)

    getGene(geneStrand,'+',upGene,geneArray, uidUP)

    return geneArray, regulatorIndex



if __name__=="__main__":
    ramr = "WP_000113609" #good 7170
    ttgr = "WP_014859138" #good 10681
    hrtr = "NP_266817.1" #good 5147
    bioq = "WP_011728885.1" #UIDs aren't neighbors
    actr = "WP_011030045.1" #good 2933
    camr = "WP_146114525.1" #No UID
    acur = "WP_011336736.1" #good 5253
    qacr = "WP_001807342.1" #UIDs aren't neighbors
    beti = "NP_414847.3" #UIDs aren't neighbors
    eilr = "WP_013366341.1" #UIDs aren't neighbors
    tetr = "WP_000113282.1" #good 3586
    bm3r1 = "WP_013083972.1" #good 4426
    pfmr = "WP_011229253.1" #UIDs aren't neighbors
    cgmr = "WP_011015249.1" #good 3353 
    cmer = "WP_002857627.1" #UIDs aren't neighbors
    sco7222 = "NP_631278.1" #good 6564
    eca1819 = "WP_011093392.1" #good 3363

    regName = bioq
    Meta = acc2MetaData(regName)
    #Meta = ['NC_003197.2', '638149', '638730', '-']
    #print(Meta)
    UID = getUID(Meta[0],Meta[1],Meta[2])
    print(UID)
    
    #linkID = getLinkID(UID)
    #print(linkID)
    
    '''
    MetaData = uid2MetaData(UID)
    
    
    time.sleep(1)
    TU, regIndex = getTU(UID, MetaData[2], MetaData[4])
    pprint(TU)
    #print(TU[regIndex])
    if abs(MetaData[2]-TU[0][2]) > 10000:
        print('something went wrong')
    else:
        print('everything checks out')

    with open('operon.data', mode='wb') as f:
        pickle.dump(TU, f)
   ''' 
