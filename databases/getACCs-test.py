import requests
import json
import time
import math
import pickle

from Bio import Entrez

Entrez.email = "doelsnitz@utexas.edu"
Entrez.api_key = "2885054140ec08a730e578c2741ff3af2008"



def searchTerms2UIDs(outfile, term):

    with open(f"{outfile}", mode="w+") as f:
        
            #count number of proteins in NCBI database that fit your query
        total_ACCs = Entrez.read(Entrez.esearch(db="protein", rettype="acc", term=term, usehistory="Y"))["Count"]
        print("total number of accessions: "+str(total_ACCs))
        
            #calculate number of requests you need to make based on above number
        numRequests = math.ceil(int(total_ACCs)/100000)

        accession_list = []

        #for i in range(0,numRequests):
        for i in range(0,1):
            time.sleep(3)
                #by default UIDs are retrieved in groups of 100K, to minimize # of requests
            retstart = i*10000

            read_success = False
            read_attempts = 0
            sleep_time = 2
            
            print("starting esearch")
            record = Entrez.read(Entrez.esearch(db="protein", retmax=1000, retstart=retstart, rettype="acc", term=term, usehistory="Y"))
            '''
            while not read_success:
                try:
                    record = Entrez.read(Entrez.esearch(db="protein", retmax=10000, retstart=retstart, rettype="acc", term=term, usehistory="Y"))
                    read_success = True
                except RuntimeError:
                    read_attempts += 1
                    print('runtime error encountered. sleeping for '+ str(sleep_time)+' seconds')
                    time.sleep(sleep_time)
                    sleep_time += 2

            if read_attempts > 10:
                print('max number or fetch attempts made. Boo!')
                break
            '''
            print("done with esearch")
            idlist = record["IdList"]

            print(idlist)
            print("starting efetch")
            h = Entrez.efetch(db="protein", id=idlist, rettype="acc")
            accString = h.read()

            #accString = "".join((item+",") for item in idlist)
        
            print("got accessions for iteration "+str(i)+" out of "+str(numRequests))

            #accession_list.append(accString)
        
        f.write(accString)
        #pickle.dump(accession_list, f)

if __name__ == "__main__":
    
    term = "txid2[ORGN] OR txid2157[ORGN] AND ( marr family OR marr repressor OR marr regulator )" 

    searchTerms2UIDs("newACCs.txt", term)
    #searchTerms2UIDs("MarR_ACCs_100K.pkl", term)


