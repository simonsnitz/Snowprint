import requests
import time
import pickle
import xmltodict

base= "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
epost = base+"epost.fcgi?db=nuccore&id="

with open("formattedTetRs.pkl", mode="rb") as f:
    uids = pickle.load(f)

response = requests.get(epost+uids[0])
if response.ok:
    data = xmltodict.parse(response.text)
    QueryKey = data["ePostResult"]["QueryKey"]
    WebEnv = data["ePostResult"]["WebEnv"]
else:
    print('bad epost request')
    print(response.status_code)


efetch = base + "efetch.fcgi?db=nuccore&rettype=acc&usehistory=y&WebEnv="+WebEnv+"&query_key="+QueryKey
reponse2 = requests.get(efetch)
if response.ok:
    print(response.text)
else:
    print('bad efetch request')
    print(response.status_code)








'''
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
'''
