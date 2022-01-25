import requests
import os
import re
import json
import time
import pickle
from Bio import Entrez

''' available functions: 

#main:
    acc2operon(protein_accession)

acc2MetaData(acc)
NC2genome(NCacc)
parseGenome(genome, seq_start, seq_stop)
fasta2MetaData(allGenes[index])
getOperon(allGenes, index, seq_start, direction)

'''    

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
        # can't use 'get' method on a list. old code = proteinList = [{}]
        proteinList = {}

    protein = proteinList.get("CDSList", "EMPTY")
    try:
        protein = protein[0].attributes
    except:
        pass

    return protein
    #[protein['accver'],protein['start'],protein['stop'],protein['strand']]


def NC2genome(NCacc):
    #this has failed in the past due to error: 'requests.exceptions.NetworkError'. Failed to establist new connection. Find fix for this?.
    
    #added sleep for 1/4 sec because I was getting 429 HTML errors that returned "efetch query unsuccessful..."
    time.sleep(0.25) 
    response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+NCacc+'&rettype=fasta_cds_aa')

    if response.ok:
        data = response.text
        with open('cache/genome.txt', mode='w+') as f:
            f.write(data)
    else:
        print(response.status_code)
        print('efetch query unsuccessful. Genome could not be found. re-trying ...')
        
        success = False
        counter = 1
        while(success == False):
            time.sleep(counter) 
            response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+NCacc+'&rettype=fasta_cds_aa')

            if response.ok:
                data = response.text
                with open('cache/genome.txt', mode='w+') as f:
                    f.write(data)
                success = True
            else:
                print('attempt number '+str(counter+1)+' failed')
                counter += 1
        

    with open('cache/genome.txt', mode='r+') as f:
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
            geneIndex += 1
            allGenes.append(i)
    try:
        return allGenes, regIndex
    except:
        print('regulator not found in genome')


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
    
    '''
    if 'link' in metaData.keys():
        try:
            #metaData['link'] = accession2UniprotID(metaData['link'])
            print('link is '+metaData['link'])
        except:
            pass
    else:
    '''
    if 'link' not in metaData.keys():
        metaData['link'] = ""   #eventually should be empty. Jinja checks if it should set a link or not
    
    return metaData

def getOperon(allGenes, index, seq_start, strand):
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
            if direction == '+':
                nextIndex = index+1
            elif direction == '-':
                nextIndex = index-1
                
            try:
                nextGene = fasta2MetaData(allGenes[nextIndex])

                if abs(seq_start - nextGene['start']) > 8000:       #added this. break if too far away
                    break
                elif geneStrand == '-' and nextGene['direction'] == '+' and direction == '+':
                    geneList.append(nextGene)
                elif geneStrand == '+' and nextGene['direction'] == '-' and direction == '-':
                    geneList.append(nextGene)
                elif geneStrand == nextGene['direction']:
                    geneList.append(nextGene)
                index = nextIndex
            except:
                break

    geneStrand = strand
    
    #attempt to get downstream genes, if there are any genes downstream
    try:
        indexDOWN = index-1
        downGene = fasta2MetaData(allGenes[indexDOWN])
        #if seq_start > downGene['start']:
        if strand == '+' and downGene['direction'] == '-':
            geneStrand = downGene['direction']
    
        downgenes = [downGene]
        getGene(geneStrand,'-',downGene, downgenes, indexDOWN)
    
        geneArray = list(reversed(downgenes))
    except:
        geneArray = []

    geneArray.append(fasta2MetaData(allGenes[index]))
    regulatorIndex = (len(geneArray)-1)

    geneStrand = strand
    
    #attempt to get upstream genes, if there are any genes upstream
    try:
        indexUP = index+1
        upGene = fasta2MetaData(allGenes[indexUP])
        #if seq_start > upGene['start']:
        if strand == '-' and upGene['direction'] == '+':
            geneStrand = upGene['direction']
        
        geneArray.append(upGene)

        getGene(geneStrand, '+', upGene, geneArray, indexUP)
    except:
        return geneArray, regulatorIndex

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

def acc2operon(accession):
    time.sleep(0.25)
    metaData = acc2MetaData(accession)
    if metaData != "EMPTY":
        time.sleep(0.25)
        genome = NC2genome(metaData["accver"])
        allGenes, index = parseGenome(genome, metaData["start"], metaData["stop"])
        reg = fasta2MetaData(allGenes[index])
        operon, regIndex = getOperon(allGenes, index, reg['start'], reg['direction'])

        data = {"operon": operon, "regIndex": regIndex, "genome": metaData["accver"], "organism":metaData["org"] }

        return data
    else:
        return "EMPTY"

#def enzyme_acc2regulator(accessions, outfile):
def enzyme_acc2regulator(accessions, max_regulators=20):
    
    regulator = re.compile(r"regulator|repressor|activator")
    operons_with_regulators = []

    if type(accessions) != list:
        accessions = accessions[:-1].split("\n")
    
    number_accessions = len(accessions)

    max_regulators = int(max_regulators)

    current_accession = 0
    total_regulators = 0

    for i in accessions:
        
        if total_regulators == max_regulators:
            return operons_with_regulators    
        
        if i[4] == "_":
            pass
        else:
            number_regulators = 0
            try:
                operon = acc2operon(i)
                #regDirection = operon["regIndex"]["direction"]
                #direction = operon["operon"][regIndex]["direction"]
                for j in operon["operon"]:
                    if regulator.search(j["description"]) and number_regulators == 0:
                    
                        #if regDirection == "-" and j["direction"] == "+":   #not working(?)
                        #    print('wrong direction')
                        #    pass
                        #else:
                        print('found regulator for '+str(i)+"! accession #: "+str(current_accession))
                        operons_with_regulators.append(operon)
                        number_regulators += 1
                        total_regulators += 1
                        print("total number of regulators = "+str(total_regulators))
                current_accession += 1
                if number_regulators == 0:
                    print('no regulator for '+str(i)+". accession #: "+str(current_accession))
            except:
                print('no operon data for '+str(i)+". accession #: "+str(current_accession))
                current_accession += 1
            

    
    return operons_with_regulators
    
    #with open(f"{outfile}", mode="wb") as out:
        #pickle.dump(operons_with_regulators, out)


if __name__ == "__main__":
    
    #enzyme_acc2regulator("WP_033004998.1")
    
    #acc2MetaData(acc)
    #NC2genome(NCacc)
    #parseGenome(genome, seq_start, seq_stop)
    #fasta2MetaData(allGenes[index])
    #getOperon(allGenes, index, seq_start, direction)
    
    #testing a BUNCH of accessions to make sure it's robust
    ramr = "WP_000113609" #good 7170                #alias good                     #genomeGood 289
    ttgr = "WP_014859138" #good 10681               #alias good                     #genomeGood 258 
    hrtr = "NP_266817.1" #good 5147                 #Not working with alias         #genomeGood 188
    bioq = "WP_011728885.1" #UIDs aren't neighbors  #alias good (fixed)             #genomeGood 111
    actr = "WP_011030045.1" #good 2933              #alias good (fixed)             #genomeGood 111
    mybCamr = "WP_146114525.1" #No UID                 #still problematic           #genomeGood 191
    acur = "WP_011336736.1" #good 5253              #RS alias good                  #genomeGood 533
    qacr = "WP_001807342.1" #UIDs aren't neighbors  #alias good                     #genomeGood 178
    beti = "NP_414847.3" #UIDs aren't neighbors     #alias good (fixed)             #genomeGood 124
    eilr = "WP_013366341.1" #UIDs aren't neighbors  #alias good (fixed)             #genomeGood 118
    tetr = "WP_000113282.1" #good 3586              #alias good                     #genomeGood 96
    bm3r1 = "WP_013083972.1" #good 4426             #alias good                     #genomeGood 304
    pfmr = "WP_011229253.1" #UIDs aren't neighbors  #alias good (fixed)             #genomeGood 308
    cgmr = "WP_011015249.1" #good 3353              #alias good                     #genomeGood 461
    cmer = "WP_002857627.1" #UIDs aren't neighbors  #alias good (fixed)             #genomeGood 159
    sco7222 = "NP_631278.1" #good 6564              #alias good                     #genomeGood 147
    eca1819 = "WP_011093392.1" #good 3363           #alias good                     #genomeGood 136
    acnr = "WP_003856101.1"                         #alias good                     #genomeGood 181
    ethr = "WP_003399797.1"  #downstream alias has "c" in name. Rv3855 -> Rv3854c   #genomeGood 76
    rutr = "WP_000191701.1"                         #alias good                     #genomeGood 231
    acrr = "WP_000101737.1"                         #alias good                     #genomeGood 142
    fadr = "WP_003229547.1" #alias number changes by 5, but doesn't have _RS        #genomeGood 134
    fadr2 = "NP_390733.1"                                                           #genomeStillGood
    marr = "WP_000843414.1"                         #alias good                     #genomeGood 212
    trpr = "WP_000068679.1"                                                         #genomeGood 211
    camr = 'BAA03510.1'                                                             #plasmidGood 323
    tcuCamr = "WP_145928353.1"                                                      #genomeGood 280
    
    '''
    data = acc2operon(camr)

    from getIntergenic import getIntergenicSeq

    intergenic = operon2Intergenic(data["operon"], data["regIndex"], data["genome"])
    print(intergenic)

    with open('cache/intergenic.txt', mode='w+') as f:
        f.write(intergenic)

    '''
    regACC = ttgr
    MetaData = acc2MetaData(regACC)
    print(MetaData["org"])
    #NCacc = MetaData[0]
    #print(NCacc)
    '''
    genome = NC2genome(MetaData[0])
    allGenes, index = parseGenome(genome, MetaData[1], MetaData[2])
    
    reg = fasta2MetaData(allGenes[index])
    operon, regIndex = getOperon(allGenes, index, reg['start'], reg['direction'])
    
    print(operon)
    if "alias" in reg.keys():
        operon[regIndex]['description'] = reg['alias'] + " regulator"

    length, intergenic = getOperator(operon, regIndex, NCacc)
    print(length, intergenic)


    with open('operon.data', mode='wb') as f:
        pickle.dump(operon, f)
    '''
