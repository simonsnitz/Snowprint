from Bio import Entrez

headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 

start = 0
finish = 3000000
esearch_inc = 100000
efetch_inc = 1000

# initialize some default parameters
Entrez.email = 'doelsnitz@utexas.edu'     # provide your email address
db = 'protein'                            # set search to dbVar database
paramEutils = { 'usehistory':'Y', 'retmax': esearch_inc }        # Use Entrez search history to cache results
term = "tetr"



with open("NCBI_TetRs.fasta", "w+") as f:

    for i in range(0, int(finish/esearch_inc)):
            # generate query to Entrez eSearch
        eSearch = Entrez.esearch(db=db, term=term, **paramEutils, retstart=i*esearch_inc)
            # get eSearch result as dict object
        res = Entrez.read(eSearch)
            # Extract IDs
        IdList = [i for i in res["IdList"]]
        num_ids = len(IdList)
        print("esearch batch "+str(i)+" of 28 complete")
        eSearch.close()

        webEnv1 = res["WebEnv"]
        queryKey1 = res["QueryKey"]

        if i > 6:

            ePost = Entrez.read(Entrez.epost(db=db, id=",".join(IdList), webenv=webEnv1, \
                query_key=queryKey1, usehistory="y"))

            webEnv2 = ePost["WebEnv"]
            queryKey2 = ePost["QueryKey"]

            for j in range(0, int(esearch_inc/efetch_inc)):

                complete = False
                while complete == False:
                    try:
                        
                        fetchHandle = Entrez.efetch(db=db, retstart=j*efetch_inc, retmax=efetch_inc, \
                            rettype="fasta", webenv=webEnv2, query_key=queryKey2, usehistory="y")
                        data=fetchHandle.read()
                        fetchHandle.close()
                        f.write(data)
                        print("efetched batch "+str(j)+" of 100 complete")
                        complete = True
                    except:
                        print("efetch failed. Retrying")
                        complete = False
                        