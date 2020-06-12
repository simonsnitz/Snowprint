import requests
import json
import time
import math
import pickle


def searchTerms2UIDs(outfile, **kwargs):

    with open(f"{outfile}", mode="w+") as f:
            
        URL_COUNT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"+\
            "db=protein&&retmode=json&rettype=count&term="
        
        organism_IDs = "txid2[ORGN]%20OR%20txid2157[ORGN]%20AND%20"    
        search_terms = "(%20marr%20family%20OR%20marr%20repressor%20OR%20marr%20regulator%20)" 
        #search_terms = "(%20tetr%20family%20OR%20tetr%20repressor%20OR%20tetr%20regulator%20)"
       
            #count number of proteins in NCBI database fit your query
        response = requests.get(URL_COUNT + organism_IDs + search_terms)

        if response.ok:
            data = response.text
            total_UIDs = json.loads(data)["esearchresult"]["count"]
            print("total number of UIDs: "+str(total_UIDs))
                #determine number of requests you'll have to make below (assuming you get 100K UIDs each request)
            numRequests = math.ceil(int(total_UIDs)/100000)
        else:
            print('Request to count number of UIDs failed')
            print(response.status_code)

        protein_list = ""

        for i in range(0,numRequests):
            time.sleep(0.5)
                
                #by default UIDs are retrieved in groups of 100K, to minimize # of requests
            retstart = str(i*100000)

            URL_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"+\
            "db=protein&retstart="+retstart+"&retmax=100000&retmode=json&term="


            response = requests.get(URL_BASE + organism_IDs + search_terms)

                #Example:
            #Archaea OR Bacteria AND (TetR family AND TetR repressor)
            #response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=txid2[ORGN]%20OR%20txid2157[ORGN]%20AND%20(%20tetr%20family%20OR%20tetr%20repressor%20OR%20tetr%20regulator%20)&retstart="+retstart+"&retmax=100000&retmode=json")

            if response.ok:
                print(str(i+1)+" out of "+str(numRequests))
                data = response.text
                idList = json.loads(data)["esearchresult"]["idlist"]
                #indexList = ""
                for idItem in idList:
                    protein_list += (idItem+"\n")
                    #f.write('%s,' % idItem)
                #protein_list.join(indexList)
            else:
                print('Request to get UIDs at iteration '+str(i)+" failed")
                print(response.status_code)
        
        f.write(protein_list)

if __name__ == "__main__":

    searchTerms2UIDs("MarR_UIDs_100K.txt")


