from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.NCBIXML import read

from pathlib import Path
from pprint import pprint

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import requests
import pickle

''' functions:
accID2sequnce - get sequence in fasta format from an input protein refseq ID (API query)
blast_and_cache - blast input protein sequence, cache output file.
homologs2accID_ident - extract accessionIDs and percent identity from homologs, store them in dictionary as output
'''


def accID2sequence(accID):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("bad request")
        print(response.status_code)


def blast_and_cache(sequence, cacheFile, perc_ident, db="nr"):
        with open(f'{cacheFile}', mode="w+") as f:
            print("found this sequence:\n"+sequence)
            print('entering blast function')
            blast_results = qblast("blastp", db, sequence, perc_ident=perc_ident)
                    #hitlist_size=alignments)
            print('finished blast')
            f.write(blast_results.read())
            print('cached blast result')
            return blast_results.read()


#def getHomologAccessionIDs(fileIN, fileOUT):
def homologs2accID_ident(fileIN, fileOUT):
    with open(f'{fileIN}', 'r') as f:
        
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

    with open(f'{fileOUT}', 'wb') as f:
        pickle.dump(homologList, f)
    
    return homologList


def acc2homolog_list(acc, blastCache, perc_ident, homologListFile):
    sequence = accID2sequence(acc)
    blast_and_cache(sequence, blastCache, perc_ident)
    homologs2accID_ident(blastCache, homologListFile)


if __name__=="__main__":

    #camr
    acc = "BAA03510.1"

    #ramr
    #acc = "WP_000113609.1"

    sequence = accID2sequence(acc)
    print(sequence)

    blast_and_cache(sequence, 'cache/blastCache/camr2000.xml')

    homologs2accID_ident('cache/blastCache/camr2000.xml', 'cache/ramr_homologList.xml')




