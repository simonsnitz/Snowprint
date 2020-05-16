import pickle
from pprint import pprint
import re
import json

with open('operon.data', 'rb') as f:
    operon = pickle.load(f)

#pprint(operon)


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
    transporter = re.compile(r"transport|pump")
    enzyme = re.compile(r"reductase|synthase|synthetase")

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
    geneEnd = max(operon[index][2],operon[index][3])
    try:
        nextGeneStart = min(operon[index+1][2],operon[index+1][3])
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
    operonStart = min(operon[0][2], operon[0][3])
    operonEnd = max(operon[-1][2], operon[-1][3])
    operonLength = operonEnd - operonStart

    if operonLength >= 6000:
        displayLength = operonLength*1.05
    else:
        displayLength = 6000
    totalLength = 0

        #populate graphic dictionary with operon gene info
    for i in range(0, len(graphic)):
        genLen = calcLength(operon[i][2], operon[i][3], displayLength)
        graphic[i]["length"] = str(genLen)+'%'
   
        graphic[i]["description"] = operon[i][0]

        geneType = type2color(operon[i][0])
        graphic[i]["color"] = geneType

        graphic[i]["direction"] = operon[i][4]

        spacer = calcSpacer(operon, i, displayLength)
        graphic[i]["spacer"] = str(spacer)+'%'
        
        #totalLength += graphic[i]["length"]
        #totalLength += graphic[i]["spacer"]

    print(operonLength)
    #print(totalLength)
    jsonGraphic = json.dumps(graphic)
    return jsonGraphic
    #return graphic

graphic = createGraphic(operon)
pprint(graphic)

with open("graphic.js", "w+") as f:
    f.write("var graphic = '"+graphic+"'")
