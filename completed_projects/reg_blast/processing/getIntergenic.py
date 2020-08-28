import pickle
import requests
from Bio import SeqIO
import processing.acc2operon as a2o

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

    if operon[regIndex]["direction"] == "-":
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
        with open('cache/intergenic.fasta', mode='w+') as f:
            f.write(intergenic)
    else:
        print('bad request')

         #800bp cutoff for an inter-operon region. A region too long makes analysis fuzzy and less accurate.
    if length <= 800:
        output = ""
        for i in intergenic.split('\n')[1:]:
            output += i
        return {"intergenicSeq": output, "regType": regType}
    else:
        print('intergenic over 800bp detected')
        return None


def appendIntergenic(homologListFile):
    with open(f'{homologListFile}', mode="rb") as f:
        homologList = pickle.load(f)

    print("original length "+str(len(homologList)))

    noDataCount = 0
    flaggedHomologs = []
    for i in range(0,len(homologList)):
        data = a2o.acc2operon(homologList[i]["accession"])
        try:
            homologList[i]["regIndex"] = data["regIndex"]
            homologList[i]["operon"] = data["operon"]
            homologList[i]["organism"] = data["organism"]           
 
            intergenic = operon2Intergenic(data["operon"], data["regIndex"], data["genome"])
            homologList[i]["intergenic"] = intergenic["intergenicSeq"]
            homologList[i]["regType"] = intergenic["regType"]
            print('got intergenic region for '+str(i)+" out of "+str(len(homologList)))
        except:
            print("no data for intergenic region at position "+str(i))
            flaggedHomologs.append(homologList[i])
            noDataCount += 1


        #create new homolog list without flagged homologs
    new_homologList = [i for i in homologList if i not in flaggedHomologs]
    print(new_homologList)

        #save updated homologListFile everytime new intergenic region appended
    with open(f'{homologListFile}', mode="wb") as f:
            pickle.dump(new_homologList,f)
            print("cached homolog file with intergenic regions appended")

    print("data not found for "+str(noDataCount)+" out of "+str(len(homologList)))




if __name__ == "__main__":

    print('main')    

