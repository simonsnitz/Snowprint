import pickle
from pprint import pprint
import re
import json

with open('operon.data', 'rb') as f:
    operon = pickle.load(f)


#calculate length of gene, place into bins
def calcLength(seq_start, seq_stop):
    size = abs(seq_start - seq_stop)
    if size > 0 and size <= 300:
        geneLen = 1
    elif  size > 500 and size <= 1000:
        genLen = 2
    elif size > 1000 and size <= 2000:
        genLen = 3
    elif size > 2000:
        genLen = 4
    return genLen

#assign colors to specific gene classes
def type2color(geneName):
    regulator = re.compile(r"regulator")
    hypothetical = re.compile(r"hypothetical")
    transporter = re.compile(r"transport")
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
        return "white"

def createGraphic(operon):
    #create empty dictionary
    graphic = {}
    for i in range(0, len(operon)):
        graphic[i] = {}

    #populate graphic dictionary with operon gene info
    for i in range(0, len(graphic)):
        genLen = calcLength(operon[i][2], operon[i][3])
        graphic[i]["length"] = genLen
   
        graphic[i]["description"] = operon[i][0]

        geneType = type2color(operon[i][0])
        graphic[i]["color"] = geneType

        graphic[i]["direction"] = operon[i][4]

    jsonGraphic = json.dumps(graphic)
    return jsonGraphic

graphic = createGraphic(operon)

with open("graphic.json", "w+") as f:
    f.write(graphic)
