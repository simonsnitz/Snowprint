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


def blast_and_cache(sequence, acc, perc_ident=50, db="nr", hitlist_size=70):
        with open(f'cache/blast_cache/{acc}.xml', mode="w+") as f:
            print("found this sequence:\n"+sequence)
            print('entering blast function')
            blast_results = qblast("blastp", db, sequence, perc_ident=perc_ident, hitlist_size=hitlist_size)
                    #hitlist_size=alignments)
            print('finished blast')
            f.write(blast_results.read())
            print('cached blast result')
            return blast_results.read()


#def getHomologAccessionIDs(fileIN, fileOUT):
def homologs2accID_ident(acc):
    with open(f'cache/blast_cache/{acc}.xml', 'r') as f:
        
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
  
    with open(f'cache/homolog_metadata/{acc}.pkl', 'wb') as f:
        pickle.dump(homologList, f)
        print('overwriting data')
    
    return homologList



def homologs2residueFrequency(acc, identity):
    with open(f'cache/blast_cache/{acc}.xml', 'r') as f:

        blast_result = read(f)
            #extract useful info from blast xml file
        metaData = [{"sequence": align.sbjct, "start": align.query_start, "end": align.query_end, "identity":round((align.identities/align.align_length)*100,2)}
                for data in blast_result.alignments
                    for align in data.hsps 
        ]
        #print(metaData)


            #create protein dataset (dictionaries within List)
        reference = metaData[0]["sequence"]
        protein = [{residue:1}
            for residue in reference
        ]

        
            #populate amino acid frequency list with homolog sequences that meet specified identity cutoff.    Maybe could use list comprehension?
        for homolog in metaData[1:]:
            if homolog["identity"] >= identity:
                counter = 0
                for pos in range(homolog["start"]-1,homolog["end"]):
                    hresidue = homolog["sequence"][counter]
                    try:
                        protein[pos][hresidue] += 1
                    except:
                        protein[pos].update({hresidue:1})
                    counter += 1

        print(protein)
        






def acc2homolog_list(acc, perc_ident):

        #check to see if blast result is already cached
    try:
        homologs = homologs2accID_ident(acc)
    except:
        print('no existing cache found')
        sequence = accID2sequence(acc)
        blast_and_cache(sequence, acc, perc_ident)
        homologs = homologs2accID_ident(acc)
    return homologs


if __name__=="__main__":

        #camr
    #acc = "BAA03510"

        #alkx
    acc = "AEM66515"

        #ramr
    #acc = "WP_000113609"

    #sequence = accID2sequence(acc)
    #print(sequence)

    homologs2residueFrequency(acc,70)

    #blast_and_cache(sequence, 'cache/blastCache/camr2000.xml')

    #homologs2accID_ident('cache/blastCache/camr2000.xml', 'cache/ramr_homologList.xml')




