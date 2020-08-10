import requests
import os

headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 


def accID2sequence(accID):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        print('good request')
        seqList = response.text.split("\n")[1:-2]
        seq =  "".join(i for i in seqList)
        return seq
    else:
        print("bad request")
        print(response.status_code)


def seq2residuePerc(seq:str,residues,name):
    seqLen = len(seq)
    count = 0
    for residue in seq:
        if residues.count(residue) > 0:
            count += 1
    percentResidue = (count/seqLen)*100
    print(name + ":   "+ str(percentResidue))


if __name__ == "__main__":
   
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


    seq = accID2sequence(bm3r1)
    aromatic = ['Y','F','W']
    hydrophobic = ['A','V','I','L','F','M','W']
    polar = ['N','D','E','R','K','S','T','H','Y','Q']
    seq2residuePerc(seq,aromatic,"aromatic")
    seq2residuePerc(seq,hydrophobic,"hydrophobic")
    seq2residuePerc(seq,polar,"polar")



