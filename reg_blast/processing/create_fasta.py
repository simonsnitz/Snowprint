import pickle
import requests
import time
import math
from os.path import exists

all_tetrs = "../cache/all_the_regulators/metadata/copy_allTetRs.pkl"
fasta_tmp = "../cache/tmp/all_rep_regs.fasta"



def accID2sequence(accID):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("bad request")
        print(response.status_code)


def check_fasta():
    with open(fasta_tmp, mode='r') as f:

        num_entries = 0
        for line in f.readlines():
            if line[0] == ">":
                num_entries += 1

        print(str(num_entries)+" entries in fasta file")
        return num_entries


def create_fasta(finish, batch_size):
    
    with open(all_tetrs, mode="rb") as d:

        data = pickle.load(d)
            
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

        # start = time.time()
        if exists(fasta_tmp):
            with open(fasta_tmp, mode="r") as f:
                seqs = f.read()
                num_entries = check_fasta()
                start_index = math.ceil(num_entries/batch_size)
        else:
            seqs = ""
            start_index = 0

        for batch in range(start_index,int(num_batches)): 
            start = time.time()
            batch = batch*batch_size
            entry_string = ",".join(EMBLs[batch:batch+batch_size])
            fastas = accID2sequence(entry_string)
            seqs += fastas
            end = time.time()
            print("appending batch "+str(int(batch/batch_size))+" out of "+str(num_batches)+" took "+str(end-start)+" seconds")
                
            with open(fasta_tmp, mode="w+") as fasta:
                fasta.write(seqs)


create_fasta("all", 300)
#check_fasta()