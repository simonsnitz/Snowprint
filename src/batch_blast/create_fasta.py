import pickle
import requests
import time
import math
from os.path import exists

# input file
all_tetrs = "../cache/all_the_regulators/metadata/copy_allTetRs.pkl"

# output file
fasta_tmp = "../cache/tmp/all_rep_regs.fasta"


# input protein accession ID, output sequence of protein in fasta format
def accID2sequence(accID):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("bad request")
        print(response.status_code)

# check how many fasta sequences are in the output file
def check_fasta(fasta_tmp):
    with open(fasta_tmp, mode='r') as f:

        num_entries = 0
        for line in f.readlines():
            if line[0] == ">":
                num_entries += 1

        print(str(num_entries)+" entries in fasta file")
        return num_entries


# create a large fasta file for making a blast database.
    # finish = index to stop at. Input a number or "all" to indicate all entries.
    # batch_size = how many protein sequences to fetch per API request (I can't do more than 300).
def create_fasta(finish, batch_size):
    
    with open(all_tetrs, mode="rb") as d:

        data = pickle.load(d)
        
        # initialize an empty array to populate with protein EMBL IDs.
        EMBLs = []
        if finish == "all":
            for entry in data:
                EMBLs.append(entry["EMBL"])

        elif type(finish) == int:        
            for entry in range(0,finish):
                EMBLs.append(data[entry]["EMBL"])
        
        else:
            raise(ValueError)

        num_batches = math.ceil(len(EMBLs)/batch_size)
        print(str(num_batches)+" batches total")

        if exists(fasta_tmp):
            with open(fasta_tmp, mode="r") as f:
                seqs = f.read()
                num_entries = check_fasta(fasta_tmp)
                start_index = math.ceil(num_entries/batch_size)
        else:
            seqs = ""
            start_index = 0

        # fetch batches of sequences (API request)
        for batch in range(start_index,int(num_batches)): 
            start = time.time()
            batch = batch*batch_size
            entry_string = ",".join(EMBLs[batch:batch+batch_size])
            fastas = accID2sequence(entry_string)
            seqs += fastas
            end = time.time()
            print("appending batch "+str(int(batch/batch_size))+" out of "+str(num_batches)+" took "+str(end-start)+" seconds")
                
            # cache result in the output file each time a batch is complete
            with open(fasta_tmp, mode="w+") as fasta:
                fasta.write(seqs)


create_fasta("all", 300)