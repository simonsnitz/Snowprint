import os
from pprint import pprint

import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO
import requests
import json
import time
import pickle
import re

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


def alias2MetaData(alias):
    
    time.sleep(0.5)
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term='+alias+'&retmode=json', headers=headers)
    if response.ok: 
        try:
            data = json.loads(response.text)["esearchresult"]["idlist"]
            uid = data[0]
        except:
            print("alias not found in esearch")
            return
            ########################################################

    else:
        print(response.status_code)
        print("alias to uid esearch failed")
        return
    
    time.sleep(0.5)
    
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id='+uid+'&retmode=json', headers=headers)
    if response.ok: 
        data = json.loads(response.text)["result"][uid]
        description = data["description"]
        gdata = data["genomicinfo"][0]

        NCacc = gdata["chraccver"]
        seq_start = gdata["chrstart"]
        seq_stop = gdata["chrstop"]
        strand = getStrand(seq_start, seq_stop)

            #get alias that works with uniprot
        alias4uniprot = alias
        re1 = re.compile(r'_RS')
        if re1.search(alias):
            newAliases = data["otheraliases"].split(',')
            for i in newAliases:
                if re1.search(i):
                    pass
                else:
                    alias4uniprot = i

        uniprotID = alias2UniprotID(alias4uniprot)
        return [description, NCacc, seq_start, seq_stop, strand, uniprotID]
    else:
        print(response.status_code)
        print("uid to meta data esummary failed")
        return


def uid2alias(uid):
    uid = str(uid)

    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id='+uid+'&retmode=json', headers=headers)
    if response.ok: 
        data = json.loads(response.text)["result"][uid]
        try:
            aliases = data["otheraliases"]
            print(aliases)
        except:
            print("elink couldn't find any accession associated with input uid")
            return
    else:
        print('bad elink request')
        return
    
    #pick alias from list. Try to avoid newer "_RS"-based aliases, since they don't work with the Uniprot API
    aliases = aliases.split(',')
    
    re1 = re.compile(r'_RS')
    alias = ""
    for i in aliases:
        if re1.search(i):
            alias = i
    if len(alias) == 0:
        alias = aliases[0]
    
    return alias
    
def alias2UniprotID(alias):

    payload = {'from': 'GENENAME',
                'to': 'ACC',
                'format': 'list',
                'query': alias,
                }

    response = requests.get('https://www.uniprot.org/uploadlists/', params=payload, headers=headers)
    if response.ok:
            #may have multiple UniprotIDs. in that case, take the first one.
        try:
            uniprotID = response.text.split('\n')[0]
        except:
            uniprotID = response.text[:-1]
        return uniprotID
    else:
        response.raise_for_status()


def getStrand(seq_start, seq_stop):
    if seq_stop > seq_start:
        return "+"
    elif seq_start > seq_stop:
        return "-"

def getOperon(alias, seq_start, strand, shiftBy):

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

        #find length of integer part in alias name
    digits = 0
    for i in reversed(alias):
        if i.isdigit():
            digits += 1
        else:
            break
    print(digits)

    def getGene(geneStrand, direction, nextGene, geneList, ALIAS, shiftBy):
        
        while geneStrand == nextGene[4]:
            
            if direction == '+':
                nextAlias = shiftAlias('+', ALIAS, shiftBy, digits)
            elif direction == '-':
                nextAlias = shiftAlias('-', ALIAS, shiftBy, digits)
                
            nextGene = alias2MetaData(nextAlias)

            if geneStrand == '-' and nextGene[4] == '+' and direction == '+':
                geneList.append(nextGene)
            elif geneStrand == nextGene[4]:
                geneList.append(nextGene)
            print(nextAlias)
            ALIAS = nextAlias
            time.sleep(0.5)
    
    def shiftAlias(direction, alias, shiftBy, digits):

        if direction == '+':
            index = int(alias[-digits:])
            index1 = str(index+shiftBy).zfill(digits)
            nextGene = alias[:-digits]+index1
        elif direction == '-':
            index = int(alias[-digits:])
            index1 = str(index-shiftBy).zfill(digits)
            nextGene = alias[:-digits]+index1
        return nextGene

    geneStrand = strand
    
    aliasDOWN = shiftAlias('-', alias, shiftBy, digits)
    downGene = alias2MetaData(aliasDOWN)
    if seq_start > downGene[2]:
        #print('downGene is behind')
        if strand == '+' and downGene[4] == '-':
            geneStrand = downGene[4]
    else:
        pass
        #print('downGene is in front')
    downgenes = [downGene]
    getGene(geneStrand,'-',downGene, downgenes, aliasDOWN, shiftBy)
    
    geneArray = list(reversed(downgenes))
    geneArray.append(alias2MetaData(alias))
    regulatorIndex = (len(geneArray)-1)

    geneStrand = strand
    
    aliasUP = shiftAlias('+', alias, shiftBy, digits)
    upGene = alias2MetaData(aliasUP)
    if seq_start > upGene[2]:
        #print('upGene is behind')
        if strand == '+' and upGene[4] == '-':
            geneStrand = upGene[4]
    else:
        pass
        #print('upGene is in front')
    geneArray.append(upGene)

    getGene(geneStrand, '+', upGene, geneArray, aliasUP, shiftBy)

    return geneArray, regulatorIndex



if __name__=="__main__":
    ramr = "WP_000113609" #good 7170                #alias good
    ttgr = "WP_014859138" #good 10681               #alias good
    hrtr = "NP_266817.1" #good 5147                 #Not working with alias - Entrez API not up to date.
    bioq = "WP_011728885.1" #UIDs aren't neighbors  #alias good (fixed)
    actr = "WP_011030045.1" #good 2933              #alias good (fixed)
    camr = "WP_146114525.1" #No UID                 #still problematic
    acur = "WP_011336736.1" #good 5253              #RS alias good
    qacr = "WP_001807342.1" #UIDs aren't neighbors  #alias good, but something is wrong
    beti = "NP_414847.3" #UIDs aren't neighbors     #alias good (fixed)
    eilr = "WP_013366341.1" #UIDs aren't neighbors  #alias good (fixed)
    tetr = "WP_000113282.1" #good 3586              #alias good
    bm3r1 = "WP_013083972.1" #good 4426             #alias good
    pfmr = "WP_011229253.1" #UIDs aren't neighbors  #alias good (fixed)
    cgmr = "WP_011015249.1" #good 3353              #alias good
    cmer = "WP_002857627.1" #UIDs aren't neighbors  #alias good (fixed)
    sco7222 = "NP_631278.1" #good 6564              #alias good
    eca1819 = "WP_011093392.1" #good 3363           #alias good

    regName = eilr
    Meta = acc2MetaData(regName)
    #Meta = ['NC_003197.2', '638149', '638730', '-']
    UID = getUID(Meta[0],Meta[1],Meta[2])
    print(UID)
    
    alias = uid2alias(UID)

    MetaData = alias2MetaData(alias)

    time.sleep(0.5)

    re1 = re.compile(r'_RS')
    if re1.search(alias):
        print('RS alias found')
        operon, regIndex = getOperon(alias, MetaData[2], MetaData[4], 5)
    else:
        print('RS alias not found')
        operon, regIndex = getOperon(alias, MetaData[2], MetaData[4], 1)

    operon[regIndex][0] = alias + " regulator"

    print(operon)
    
    if abs(MetaData[2]-operon[0][2]) > 10000:
        print('something went wrong')
    else:
        print('everything checks out')

    with open('operon.data', mode='wb') as f:
        pickle.dump(operon, f)
   

    '''
    MetaData = uid2MetaData(UID)
    
    
    time.sleep(0.5)
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
