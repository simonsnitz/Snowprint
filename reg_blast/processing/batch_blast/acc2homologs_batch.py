from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.NCBIXML import read, parse

import requests
import pickle
import time
import os

''' functions:
accID2sequnce - get sequence in fasta format from an input protein refseq ID (API query)
create_fasta - create a fasta file in /tmp with X number of sequences from the regulator database
blast_fasta - blast input fasta file, cache output file.
append_blast_fasta - add cached blast data to the corresponding regulator in the database
'''

all_tetrs = "../cache/all_the_regulators/metadata/filtered_TetRs.pkl"
fasta_tmp = "../cache/tmp/regulators_tmp.fasta"
blast_cache = "../cache/tmp/blast_tmp.xml"




def accID2sequence(accID):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("bad request")
        print(response.status_code)




def create_fasta(finish, step):
    
    with open(all_tetrs, mode="rb") as d:

        data = pickle.load(d)
            
        EMBLs = []
        indices = []
        
        c = 0
        for entry in range(0,finish):
            if c < step and "homologs" not in data[entry].keys():
                EMBLs.append(data[entry]["EMBL"])
                indices.append(entry)
                c += 1


        with open(fasta_tmp, mode="w+") as fasta:

            seqs = [ accID2sequence(i) for i in EMBLs]
            seqs = [i for i in seqs if type(i) == str]

            s = "".join(seqs)
            print("cached fasta file")
            fasta.write(s)

    return indices




def blast_fasta():
    
    with open('../cache/tmp/blast_tmp.xml', mode="w+") as f:
        print('entering blast function')
            
        fasta = open(fasta_tmp).read()

        blastStart = time.time()

        blast_results = qblast("blastp", "nr", fasta)
            
        blastEnd = time.time()

        print('finished blast. Took '+str(blastEnd - blastStart)+" seconds.")
        f.write(blast_results.read())
        print('cached blast result')   




def append_blast_data(indices):

    with open(all_tetrs, mode="rb") as f:
        database = pickle.load(f)

    with open(blast_cache, mode="r") as f:

        records = [record for record in parse(f)]

        if len(records) > 0:
            print("there are this many blast records: "+str(len(records)))
            c = 0
            for i in indices:
                record = records[c]
                homologs = [{"accession":alignment.accession, "identity":round((hsp.identities/hsp.align_length)*100, 2)}
                    for alignment in record.alignments
                        for hsp in alignment.hsps
                ]
                homologs = [homolog for homolog in homologs if homolog["identity"] >= 50]
                
                print("the length of homolog array: "+str(len(homologs)))
                if "homologs" not in database[i].keys() and len(homologs) >= 1:
                    database[i]["homologs"] = homologs
                    print("appended homolog data for entry "+str(i))
                else:
                    print("blast did not return any alignments")
                c += 1

                #update database
            with open(all_tetrs, mode="wb") as db:
                database = pickle.dump(database, db)   
        else:
            print("there was an issue with blast") 




def check_blast():
        with open(blast_cache, mode="r") as f:

            records = [record for record in parse(f)]
            for record in records:
                homologs = [{"accession":alignment.accession, "identity":round((hsp.identities/hsp.align_length)*100, 2)}
                    for alignment in record.alignments
                        for hsp in alignment.hsps
                ]
            print(homologs)



    # compiler function
def acc2homologs(end, step):
    with open(all_tetrs, mode="rb") as db:
        database = pickle.load(db) 
    
    while "homologs" not in database[end].keys():
        
        blast_success = False

        indices = create_fasta(end, step)
        print("blasting indices: "+str(indices))
        blast_fasta()

        while blast_success == False:

            filesize = os.path.getsize(blast_cache)
            print('blast cache file size: '+str(filesize))

            if filesize > 6000:
                append_blast_data(indices)
                blast_success = True
                
            else:
                step = step-1
                print('blast failed. Re-trying with step = '+str(step))
                indices = create_fasta(end, step)
                print("blasting indices: "+str(indices))
                blast_fasta()
    



if __name__=="__main__":

    end = 1000   # what index to stop blasting at
    step = 4    # how many sequences to blast at a time

    acc2homologs(end, step)


    # should add a function to reduce the 'step' value if the blast fails