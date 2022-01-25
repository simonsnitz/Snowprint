import pickle
import requests
import src.acc2operon_batch as a2os
import time

import os


#TODO
#When running program for AlkU, which has an operator between same-direction genes, this program did not extract the correct inter-operon region. This region is ~310bp, which is well above my supposed "100bp" cutoff to find interoperon regions between same direction genes.
#Two possible causes:
#One: KEGG annotations are different than NCBI annotations. NCBI may indicate that there is a gene in that region
#Two: My program isn't catching it. Not using this filter for some reason, in this particular case.

#def getIntergenicSeq(operon, regIndex, NCacc):
def operon2Intergenic(operon, regIndex, NCacc):

    if operon[regIndex]["direction"] == "+":
        queryGenes = reversed(operon[0:regIndex])
        index = regIndex
        for i in queryGenes:
            if i["direction"] == "-":
                startPos = i["stop"]
                stopPos = operon[index]["start"]
                regType = "type 1"
                break
            else:
                start = operon[regIndex-1]["stop"]
                stop = operon[regIndex]["start"]
                testLength = int(stop) - int(start)
                    #setting this to 100bp is somewhat arbitrary. Most intergenic regions >= 100bp. May need to tweak.
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    regType = "type 2"
                    print('operator in between same-direction genes')
                    break
                index -= 1

    elif operon[regIndex]["direction"] == "-":
        queryGenes = operon[regIndex+1:]
        index = regIndex
        for i in queryGenes:
            if i["direction"] == "+":
                stopPos = i["start"]
                startPos = operon[index]["stop"]
                regType = "type 1"
                break
            else:
                    #counterintunitive use of "stop"/"start", since start < stop always true, regardless of direction
                start = operon[regIndex]["stop"]
                stop = operon[regIndex+1]["start"]
                testLength = int(stop) - int(start)
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    regType = "type 2"
                    print('operator in between same-direction genes')
                    break
                index += 1

    length = int(stopPos) - int(startPos)
  
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="+str(NCacc)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos)+"&strand=1&rettype=fasta"
    response = requests.get(URL)

    if response.ok:
        intergenic = response.text
    else:
        print('bad request')

         #800bp cutoff for an inter-operon region. A region too long makes analysis fuzzy and less accurate.
    output  = ""
    for i in intergenic.split('\n')[1:]:
        output += i
    if len(output) <= 800:
        return {"intergenicSeq": output, "regType": regType}
    else:
        print('intergenic region over 800bp')
        return None



def appendIntergenic_batch(acc):

    inputPath = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', f"cache/homolog_metadata/{acc}.pkl"))
    updatedPath = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', f"cache/homolog_metadata/updated_metadata/{acc}.pkl"))


    try:
        with open(updatedPath, "rb") as f:
            homologList = pickle.load(f)
    except:
        with open(inputPath, "rb") as f:
            homologList = pickle.load(f)


    operonProcessingTimes = [] 

    if "regType" not in homologList[0]:

        getOperons_start = time.time()

            # get operon

        accessions = [i["accession"] for i in homologList]

        data = a2os.batch_acc2operon(accessions)
        
        for i in range(0, len(data)):
            try:
                homologList[i]["regIndex"] = data[i]["data"]["regIndex"]
                homologList[i]["operon"] = data[i]["data"]["operon"]           
                        # efetch request
                intergenic = operon2Intergenic(data[i]["data"]["operon"],\
                     data[i]["data"]["regIndex"], data[i]["data"]["genome"])
                homologList[i]["intergenic"] = intergenic["intergenicSeq"]
                homologList[i]["regType"] = intergenic["regType"]
                    
                print('got intergenic region for '+str(i))

                    #cache results
                with open(updatedPath, "wb") as f:
                    pickle.dump(homologList,f)
                
            except:
                print("no data for intergenic region at position "+str(i))

        getOperons_end = time.time()

        print("operon "+str(i)+" took "+str(getOperons_end - getOperons_start)+" seconds")
        operonProcessingTimes.append(getOperons_end-getOperons_start)
        

            # remove entries without a set intergenic region
        new_homologList = [i for i in homologList if "intergenic" in i]
        with open(updatedPath, "wb") as f:
            pickle.dump(new_homologList,f)

            # print statistics
        print("total time spent getting operons: "+str(sum(operonProcessingTimes)))
        #try:
        #    print("average processing time: "+str(statistics.mean(operonProcessingTimes)))
        #    print("median processing time: "+str(statistics.median(operonProcessingTimes)))
        #except:
        #    pass

    
    else:
        print("intergenic data already cached for "+str(acc))

    


if __name__ == "__main__":

    appendIntergenic_batch("WP_000113282")

