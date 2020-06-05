from Bio import Entrez
import pickle
import time
import sys

Entrez.email = "doelsnitz@utexas.edu"

def UIDs2fastaDB(infile, outfile):

    with open(f"{infile}", mode="rb") as infile:
        UIDs = pickle.load(infile)

    numRequests = len(UIDs)


    with open(f"{outfile}", mode="w+") as outfile:

        for i in range(0, numRequests):
            time.sleep(0.5)

            post = Entrez.read(Entrez.epost(db="protein", id=UIDs[i]))
            webenv = post["WebEnv"]
            query_key = post["QueryKey"]

            print("Posted 100K IDs for iteration "+str(i+1)+" out of "+str(numRequests))
            
            time.sleep(0.5)
            accessions = Entrez.efetch(webenv=webenv,
                                        query_key=query_key,
                                        db="protein",
                                        rettype="fasta",
                                        retmode="text",
                                        usehistory='Y').read()

            print('Got accessions for iteration '+str(i+1))

            outfile.write(accessions)

if __name__ == "__main__":

    #UIDs2fastaDB("MarR_UIDs_100K.pkl", "MarR_Database.fsa")

    infile = sys.argv[1]

    outfile = sys.argv[2]

    UIDs2fastaDB(infile, outfile)
