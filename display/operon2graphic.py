import pickle
from pprint import pprint
import re
import json


    #calculate length of gene, return % size compared to total operon length
def calcLength(seq_start, seq_stop, operonLength):
    size = abs(seq_start - seq_stop)
    geneLen = round((size/operonLength)*100,0)
    geneLen = geneLen
    return geneLen


    #assign colors to specific gene classes
def type2color(geneName):
    regulator = re.compile(r"regulator|repressor")
    hypothetical = re.compile(r"hypothetical")
    transporter = re.compile(r"port|pump")
    enzyme = re.compile(r"ase|P450")

    if regulator.search(geneName):
        #regulator = green
        return "#6eff42"
    elif hypothetical.search(geneName):
        #hypothetical = grey
        return "#bfbfbf"
    elif transporter.search(geneName):
        #transporter = yellow
        return "#ffe042"
    elif enzyme.search(geneName):
        #enzyme = red
        return "#ff4242"
    else:
        #other = white
        return "black"


    #calculate the spacer length between a gene and its downstream neighbor
def calcSpacer(operon, index, operonLength):
    geneEnd = max(operon[index]['start'],operon[index]['stop'])
    try:
        nextGeneStart = min(operon[index+1]['start'],operon[index+1]['stop'])
    except:
        nextGeneStart = geneEnd
    spacer = nextGeneStart - geneEnd
    if spacer < 0:
        spacer = 0
    spacerSize = round((spacer/operonLength)*100,0)
    spacerSize = spacerSize
    return spacerSize



def createGraphic(operon):
        #create empty dictionary
    graphic = {}
    for i in range(0, len(operon)):
        graphic[i] = {}
   

        #calculate total length of operon
    operonStart = min(operon[0]['start'], operon[0]['stop'])
    operonEnd = max(operon[-1]['start'], operon[-1]['stop'])
    operonLength = operonEnd - operonStart

    if operonLength >= 6000:
        displayLength = operonLength*1.05
    else:
        displayLength = 6500
    totalLength = 0

        #populate graphic dictionary with operon gene info
    for i in range(0, len(graphic)):
        genLen = calcLength(operon[i]['start'], operon[i]['stop'], displayLength)
        graphic[i]["length"] = str(genLen)+'%'
   
        graphic[i]["description"] = operon[i]['description']

        geneType = type2color(operon[i]['description'])
        graphic[i]["color"] = geneType

        graphic[i]["direction"] = operon[i]['direction']

        spacer = calcSpacer(operon, i, displayLength)
        graphic[i]["spacer"] = str(spacer)+'%'
        
        graphic[i]["link"] = operon[i]['link']

        #totalLength += graphic[i]["length"]
        #totalLength += graphic[i]["spacer"]

    #print(operonLength)
    #print(totalLength)
    jsonGraphic = json.dumps(graphic)
    return jsonGraphic

#def saveGraphic(infile, outfile):
def saveGraphic(operons, outfile):
    
    #with open(f'{infile}', 'rb') as f:
    #operons = pickle.load(f)

    graphic = [ createGraphic(i["operon"]) for i in operons ]
    
    #pprint(graphic[1])

    with open(f"{outfile}", "w+") as f:
        f.write("var graphic = ["+str(graphic)[1:-1]+"]")


if __name__ == "__main__":

    inputFile = "octanol_operons2.pkl"
    outputFile = "octanol_graphic_native.js"

    saveGraphic(inputFile, outputFile)


