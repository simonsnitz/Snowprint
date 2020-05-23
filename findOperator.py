import pickle
import requests


def getOperator(operon, regIndex, NCacc):

    if operon[regIndex]["direction"] == "+":
        queryGenes = reversed(operon[0:regIndex])
        index = regIndex
        for i in queryGenes:
            if i["direction"] == "-":
                startPos = i["stop"]
                stopPos = operon[index]["start"]
                break
            else:
                start = operon[regIndex-1]["stop"]
                stop = operon[regIndex]["start"]
                testLength = int(stop) - int(start)
                    #setting this to 100bp is somewhat arbitrary. Most intergenic regions >= 100bp. May need to tweak.
                if testLength > 100:
                    startPos = start
                    stopPos = stop
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
                break
            else:
                    #counterintunitive use of "stop"/"start", since start < stop always true, regardless of direction
                start = operon[regIndex]["stop"]
                stop = operon[regIndex+1]["start"]
                testLength = int(stop) - int(start)
                if testLength > 100:
                    startPos = start
                    stopPos = stop
                    print('operator in between same-direction genes')
                    break
                index += 1

    length = int(stopPos) - int(startPos)
   
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="+str(NCacc)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos)+"&strand=1&rettype=fasta"
    response = requests.get(URL)

    if response.ok:
        intergenic = response.text
        with open('intergenic.fasta', mode='w+') as f:
            f.write(intergenic)
    else:
        print('bad request')

    return length, intergenic


if __name__ == "__main__":

    with open('operon.data', 'rb') as f:
        data = pickle.load(f)

    print(data)

    regIndex = 2
    operatorRegion = getOperator(data, regIndex)

    print(operatorRegion)