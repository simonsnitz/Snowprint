from Bio import Entrez
import pickle
import time
import sys

Entrez.email = "doelsnitz@utexas.edu"

def UIDs2fastaDB(infile, outfile):

    with open(f"{infile}", mode="rb") as infile:
        UIDs = pickle.load(infile)

    print(len(UIDs[0].split(',')))
    numRequests = len(UIDs)


    with open(f"{outfile}", mode="w+") as outfile:

        allAccessions = ""
        for i in range(0, numRequests):
            time.sleep(0.5)

            post = Entrez.read(Entrez.epost(db="protein", id=UIDs[i], usehistory="y"))
            webenv = post["WebEnv"]
            query_key = post["QueryKey"]

            print("Posted 100K IDs for iteration "+str(i+1)+" out of "+str(numRequests))
            
            time.sleep(1)
            for j in range(0,10):
                print("getting accessions now")
                accessions = Entrez.efetch(webenv=webenv,
                                        query_key=query_key,
                                        db="protein",
                                        rettype="fasta",
                                        retstart = (j*10000),
                                        retmax=10000,
                                        retmode="text",
                                        usehistory='y').read()
                print('done getting accessions')

                print('Got accessions for sub-iteration '+str(j+1)+"/10 of main iteration "+str(i+1)+"/"+str(numRequests))


                analyze = accessions.split('\n')
                counter = 0
                for i in analyze:
                    try:
                        if i[0] == ">":
                            counter += 1
                    except: 
                        pass
                print(counter)

                allAccessions += accessions

                outfile.write(allAccessions)

if __name__ == "__main__":

    #UIDs2fastaDB("MarR_UIDs_100K.pkl", "MarR_Database.fsa")

    infile = sys.argv[1]

    outfile = sys.argv[2]

    UIDs2fastaDB(infile, outfile)
