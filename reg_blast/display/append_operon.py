import re
import json
from pathlib import Path

    # Calculate length of gene, return % size compared to total operon length
def calcLength(seq_start, seq_stop, operonLength):
    size = abs(seq_start - seq_stop)
    geneLen = round((size/operonLength)*100,0)
    geneLen = geneLen
    return geneLen


    # Assign colors to specific gene classes
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


    # Calculate the spacer length between a gene and its downstream neighbor
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



def createGraphic(operon, targetIndex):
        
        # Create empty dictionary
    graphic = {}
    for i in range(0, len(operon)):
        graphic[i] = {}
   

        # Calculate total length of operon
    operonStart = min(operon[0]['start'], operon[0]['stop'])
    operonEnd = max(operon[-1]['start'], operon[-1]['stop'])
    operonLength = operonEnd - operonStart

    if operonLength >= 6000:
        displayLength = operonLength*1.05
    else:
        displayLength = 6500

        # Populate graphic dictionary with operon gene info
    for i in range(0, len(graphic)):
        genLen = calcLength(operon[i]['start'], operon[i]['stop'], displayLength)
        graphic[i]["length"] = str(genLen)+'%'
   
        graphic[i]["description"] = operon[i]['description']

        geneType = type2color(operon[i]['description'])
        
        if i == targetIndex:
            graphic[i]["color"] = "blue"
        else:
            graphic[i]["color"] = geneType

        graphic[i]["direction"] = operon[i]['direction']

        spacer = calcSpacer(operon, i, displayLength)
        graphic[i]["spacer"] = str(spacer)+'%'
        
        graphic[i]["link"] = operon[i]['link']

    jsonGraphic = json.dumps(graphic)
    return jsonGraphic


def create_operon_html(operons):
    

    graphic = [{"meta_data": json.dumps({"accession":i["accession"],"identity":i["identity"], \
        "organism":i["organism"]}), "operon_data": createGraphic(i["operon"], i["regIndex"])} \
            for i in operons ]


    graphic_js_variable = "var operon_graphic = ["+str(graphic)[1:-1]+"]"


    JS_string =   """

            // Operon object //

            """+graphic_js_variable+"""

            var len = operon_graphic.length
            for (i in range(0, len)){
                var metaData = JSON.parse(operon_graphic[i]["meta_data"]);
                var operon = JSON.parse(operon_graphic[i]["operon_data"]);
                renderOperon(metaData, operon);
            }
                    """

        # TODO: CHANGE THIS
    p = Path("/home/simonsnitz/projects/tuSeek/reg_blast/display/Webpage")
    fp = p / "operon_data.js"


       # Overwrite javascript file.
    with open(str(fp), mode="w") as f:
        f.write(JS_string)


