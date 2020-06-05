import requests
import json
import time
import pickle

with open("formattedTetRs.pkl", mode="rb") as f:
    accList = pickle.load(f)


with open("TetRs_extras.fsa", mode="w+") as f:
    #for i in range(7455, len(accList)):
    for i in [2653, 2836, 2837, 3340, 3345, 3349, 3354, 3348, 3355, 3361, 3364, 3365, 3367, 5311, 6247, 6795, 7086, 7216, 7997]:
        time.sleep(0.1)

        #Archaea OR Bacteria AND (TetR family AND TetR repressor)
        try:
            response = requests.post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+accList[i]+"&rettype=fasta&retmode=text")

            if response.ok:
                try:
                    data = response.text
                    f.write(data)
                    print("saved "+str(i)+" of 7997")
                except:
                    print("there was an error at iteration "+str(i))
            else:
                print(response.status_code)
                print("there was an error at iteration "+str(i))
        except:
            print('there was a connection error at iteration '+str(i))
            continue
