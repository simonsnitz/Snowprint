from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.NCBIXML import read

from pathlib import Path
from pprint import pprint

#simon added this
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import requests
import pickle

def accID2sequence(accID):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("bad request")
        print(response.status_code)


def blast_and_cache(sequence, db="nr", alignments=50):
        with open(f'results.xml', mode="w+") as f:
            print('entering blast function')
            blast_results = qblast("blastp", db, sequence, hitlist_size=alignments)
            print('finished blast')
            f.write(blast_results.read())
            print('cached blast result')


def openFasta(filePath: str) -> SeqRecord:
    with open(filePath, mode = "r") as f:
        return SeqIO.read(f, 'fasta')


def getHomologAccessionIDs(filename):
    with open(f'{filename}', 'r') as f:
        
        blast_results = read(f)
        '''     more readable, but slower
        def calcIdentity(alignment):
            length = alignment.hsps[0].align_length
            matches = alignment.hsps[0].identities
            identity = round((matches/length)*100,2)
            return identity

        homologList = []

        for alignment in blast_results.alignments:
            accession = alignment.accession
            identity = calcIdentity(alignment)
            homologList.append({"accession":accession,"identity":identity})
        '''
        #less readable but faster
    homologList = [{"accession":alignment.accession, "identity":round((hsp.identities/hsp.align_length)*100, 2)}
            for alignment in blast_results.alignments
                for hsp in alignment.hsps
         ]
  
    #reduce homolog list (so you don't have to wait forever)
    #homologList = homologList[0:10]

    with open(f'homolog_list.pkl', 'wb') as f:
        pickle.dump(homologList, f)
    
    return homologList


if __name__=="__main__":


    acc = "BAA03510.1"

    #sequence = accID2sequence(acc)

    #blast_and_cache(sequence)

    pprint(getHomologAccessionIDs('results.xml'))


    '''
    E_VALUE_THRESH = 0.1
    for alignment in blast_results.alignments:
        print(alignment.__dict__)
        for hsp in alignment.hsps:
            print(hsp.__dict__)
            if hsp.expect < E_VALUE_THRESH:
                print("****Alignment****")
                print("Description: ", alignment.title)
                print("Accesion: ", alignment.accession)
                print("length:", alignment.length)
                print("e value:", hsp.expect)
                print('score : ', hsp.score)
                print('Sequence similarity: ',round(hsp.identities/alignment.length, 3))
                # print(hsp.query[0:75] + "...")
                # print(hsp.match[0:75] + "...")
                # print(hsp.sbjct[0:75] + "...")
    '''



