import os
from pprint import pprint

import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO
import requests
import json
import time

Entrez.email = os.getenv("EMAIL", "danny.diaz@utexas.edu")
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
        return [description, NCacc, seq_start, seq_stop, strand]
    else:
        print(response.status_code)


def getStrand(seq_start, seq_stop):
    if seq_stop > seq_start:
        return "+"
    elif seq_start > seq_stop:
        return "-"

def getTU(uid, seq_start, strand):
    
    def getGene(geneStrand, direction, nextGene, geneList, UID):
        print('starting '+direction)
        while geneStrand == nextGene[4]:
            print(geneStrand, nextGene[4])
            if direction == '+':
                UID += 1
            elif direction == '-':
                UID -= 1
            nextGene = uid2MetaData(UID)
            if geneStrand == nextGene[4]:
                geneList.append(nextGene)
            time.sleep(1)
    
    uidUP = uid+1
    upGene = uid2MetaData(uidUP)
    geneStrand = strand
    if seq_start > upGene[2]:
        print('upGene is behind')
        if strand == '+' and upGene[2] == '-':
            geneStrand = upGene[4]
    else:
        print('upGene is in front')
    upgenes = [upGene]
    getGene(geneStrand,'+',upGene,upgenes, uidUP)

    uidDOWN = uid-1
    downGene = uid2MetaData(uidDOWN)
    if seq_start > downGene[2]:
        print('downGene is behind')
        if strand == '+' and upGene[2] == '-':
            geneStrand = upGene[4]
    else:
        print('downGene is in front')
    downgenes = [downGene]
    getGene(geneStrand,'-',downGene, downgenes, uidDOWN)

    return upgenes, downgenes



if __name__=="__main__":
    ramr = "WP_000113609"
    ttgr = "WP_014859138"
    hrtr = "NP_266817.1"
    bioq = "WP_011728885.1"
    actr = "WP_011030045.1"
    camr = "WP_146114525.1"
    acur = "WP_011336736.1"


    Meta = acc2MetaData(acur)
    #Meta = ['NC_003197.2', '638149', '638730', '-']
    print(Meta)
    UID = getUID(Meta[0],Meta[1],Meta[2])
    print(UID)
    MetaData = uid2MetaData(UID)
    #print(MetaData)
    
    time.sleep(1)
    TU = getTU(UID, MetaData[2], MetaData[4])
    print(TU)
    if abs(MetaData[2]-TU[0][0][2]) > 10000:
        print('something went wrong')
    else:
        print('everything checks out')

