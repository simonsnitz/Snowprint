"""
***ABOUT***

    This is a prototype version of acc2operon.

    An improved verion of these functions is found in the 'acc2operon.py' file.
"""

import requests
import os
import re
import json
import time
import pickle
from Bio import Entrez

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
    print(protein)

    return [protein['accver'],protein['start'],protein['stop'],protein['strand']]


def NC2genome(NCacc):
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+NCacc+'&rettype=fasta_cds_aa')

    if response.ok:
        data = response.text
        with open('genome.txt', mode='w+') as f:
            f.write(data)
            print('got genome')
    else:
        print('bad reqiuest')
    
    with open('genome.txt', mode='r+') as f:
        genome = f.readlines()
    return genome


def parseGenome(genome, start, stop):
    re1 = re.compile(start)
    re2 = re.compile(stop)
    geneIndex = 0
    allGenes = []
    for i in genome:
        if i[0] == '>':
            if re1.search(i):
                if re2.search(i):
                    regIndex = geneIndex
                    print(i)
                    print('regulator found at index '+str(regIndex))
            geneIndex += 1
            allGenes.append(i)
    return allGenes, regIndex


def fasta2MetaData(fasta):
    metaData = {}
    regulator = fasta.split(' [')
    
    for i in regulator:
        if i[:10] == 'locus_tag=':
            metaData['alias'] = i[10:-1]
        elif i[:8] == 'protein=':
            metaData['description'] = i[8:-1].replace("'", "")
        elif i[:11] == 'protein_id=':
            metaData['link'] = i[11:-1]
            print(metaData['link'])
        elif i[:9] == 'location=':
            if i[9:20] == 'complement(':
                metaData['direction'] = '-'
                location = i[20:-2]
                location = location.split('..')
                metaData['start'] = int(re.sub("\D", "", location[0]))
                metaData['stop'] = int(re.sub("\D", "", location[1]))
            else:
                metaData['direction'] = '+'
                location = i[9:-1]
                location = location.split('..')
                metaData['start'] = int(re.sub("\D", "", location[0]))
                metaData['stop'] = int(re.sub("\D", "", location[1]))
    
    print('getting uniprotID')
    if 'link' in metaData.keys():
        try:
            #metaData['link'] = accession2UniprotID(metaData['link'])
            print('link is '+metaData['link'])
        except:
            pass
    else:
        metaData['link'] = ""   #eventually should be empty. Jinja checks if it should set a link or not
    
    return metaData

def getOperon(allGenes, index, seq_start, strand):
    #alias -> index
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

    def getGene(geneStrand, direction, nextGene, geneList, index):
        
        while geneStrand == nextGene['direction']:
            #nextAlias = nextIndex
            if direction == '+':
                nextIndex = index+1
                #shiftAlias('+', ALIAS, shiftBy, digits)
            elif direction == '-':
                nextIndex = index-1
                #shiftAlias('-', ALIAS, shiftBy, digits)
                
            nextGene = fasta2MetaData(allGenes[nextIndex])

            if geneStrand == '-' and nextGene['direction'] == '+' and direction == '+':
                geneList.append(nextGene)
            elif geneStrand == nextGene['direction']:
                geneList.append(nextGene)
            print(nextGene['alias'])
            index = nextIndex
            time.sleep(0.1)

    geneStrand = strand
    
        #aliasDOWN = indexDOWN
    indexDOWN = index-1
    downGene = fasta2MetaData(allGenes[indexDOWN])
        ################
    if seq_start > downGene['start']:
        #print('downGene is behind')
        if strand == '+' and downGene['direction'] == '-':
            geneStrand = downGene['direction']
    else:
        pass
        #print('downGene is in front')
    downgenes = [downGene]
    getGene(geneStrand,'-',downGene, downgenes, indexDOWN)
    
    geneArray = list(reversed(downgenes))
    geneArray.append(fasta2MetaData(allGenes[index]))
    regulatorIndex = (len(geneArray)-1)

    geneStrand = strand
    
    indexUP = index+1
    upGene = fasta2MetaData(allGenes[indexUP])
    if seq_start > upGene['start']:
        #print('upGene is behind')
        if strand == '+' and upGene['direction'] == '-':
            geneStrand = upGene['direction']
    else:
        pass
        #print('upGene is in front')
    geneArray.append(upGene)

    getGene(geneStrand, '+', upGene, geneArray, indexUP)

    return geneArray, regulatorIndex


def accession2UniprotID(alias):

    payload = {'from': 'P_REFSEQ_AC',
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


if __name__ == "__main__":
    
    #acc2MetaData(acc)
    #NC2genome(NCacc)
    #parseGenome(genome, regulator)
    #fasta2MetaData(allGenes[index])
    #getOperon(allGenes, index, seq_start, strand)
    
    ramr = "WP_000113609" #good 7170                #alias good                     #genomeGood
    ttgr = "WP_014859138" #good 10681               #alias good                     #genomeGood  
    hrtr = "NP_266817.1" #good 5147                 #Not working with alias         #genomeGood
    bioq = "WP_011728885.1" #UIDs aren't neighbors  #alias good (fixed)             #genomeGood, but missing genes
    actr = "WP_011030045.1" #good 2933              #alias good (fixed)             #genomeGood, but missing genes
    camr = "WP_146114525.1" #No UID                 #still problematic              #genomeGood
    acur = "WP_011336736.1" #good 5253              #RS alias good                  #genomeGood
    qacr = "WP_001807342.1" #UIDs aren't neighbors  #alias good                     #genomeGood
    beti = "NP_414847.3" #UIDs aren't neighbors     #alias good (fixed)             #genomeGood
    eilr = "WP_013366341.1" #UIDs aren't neighbors  #alias good (fixed)             #genomeGood
    tetr = "WP_000113282.1" #good 3586              #alias good                     #genomeGood
    bm3r1 = "WP_013083972.1" #good 4426             #alias good                     #genomeGood
    pfmr = "WP_011229253.1" #UIDs aren't neighbors  #alias good (fixed)             #genomeGood
    cgmr = "WP_011015249.1" #good 3353              #alias good                     #genomeGood, but missing gene(s). Edit operon finder
    cmer = "WP_002857627.1" #UIDs aren't neighbors  #alias good (fixed)             #genomeGood
    sco7222 = "NP_631278.1" #good 6564              #alias good                     #genomeGood
    eca1819 = "WP_011093392.1" #good 3363           #alias good                     #genomeGood
    acnr = "WP_003856101.1"                         #alias good                     #genomeGood, but missing gene(s)
    ethr = "WP_003399797.1"  #downstream alias has "c" in name. Rv3855 -> Rv3854c   #genomeGood
    rutr = "WP_000191701.1"                         #alias good                     #genomeGood
    acrr = "WP_000101737.1"                         #alias good                     #genomeGood
    fadr = "WP_003229547.1" #alias number changes by 5, but doesn't have _RS        #genomeGood
    fadr2 = "NP_390733.1"                                                           #genomeStillGood
    marr = "WP_000843414.1"                         #alias good                     #genomeGood
    trpr = "WP_000068679.1"                                                         #genomeGood, but missing gene(s)


    regACC = ramr
    MetaData = acc2MetaData(regACC)
    print(MetaData[0])
    genome = NC2genome(MetaData[0])
    allGenes, index = parseGenome(genome, MetaData[1], MetaData[2])
    
    reg = fasta2MetaData(allGenes[index])
    operon, regIndex = getOperon(allGenes, index, reg['start'], reg['direction'])
    
    print(operon)
    #operon[regIndex]['alias'] = 

    with open('operon.data', mode='wb') as f:
        pickle.dump(operon, f)
