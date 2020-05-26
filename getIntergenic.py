import pickle
import requests
from Bio import SeqIO
import acc2operon as a2o


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

    if length <= 500:
        output = ""
        for i in intergenic.split('\n')[1:]:
            output += i
        return {"intergenicSeq": output, "regType": regType}
    else:
        print('intergenic over 500bp detected')
        return None


def appendIntergenic(homologListFile):
    with open(f'{homologListFile}', mode="rb") as f:
        homologList = pickle.load(f)

    for i in range(0,len(homologList)):
        protein = homologList[i]
        data = a2o.acc2operon(protein["accession"])
        try:
            intergenic = operon2Intergenic(data["operon"], data["regIndex"], data["genome"])
            homologList[i]["intergenic"] = intergenic["intergenicSeq"]
            homologList[i]["regType"] = intergenic["regType"]
        except:
            print("no data for intergenic region")
        print('got intergenic region')
        print(intergenic)
            #save updated homologListFile everytime new intergenic region appended
        with open(f'{homologListFile}', mode="wb") as f:
            pickle.dump(homologList,f)


if __name__ == "__main__":

    print('main')    

