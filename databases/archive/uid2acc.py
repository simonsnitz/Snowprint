import requests
import time
import pickle

URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=acc&id="

with open("formattedTetRs.pkl", mode="rb") as f:
    uids = pickle.load(f)

uids2 = uids[2168:]

with open("tetrACCs_2168.txt", mode="w+") as outfile:
    counter = 2168
    for line in uids2:
        time.sleep(0.1)
        response = requests.get(URL+line)
        if response.ok:
            outfile.write(response.text[:-1])
            print("got "+str(counter)+" out of "+str(len(uids)))

        else:
            print(response.status_code)
            print('bad request at iteration '+str(counter))
        counter += 1
