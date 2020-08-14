from acc2homologs import acc2homolog_list
from getIntergenic import appendIntergenic
from operator_utils import appendOperatorMetadata
from display.operator_graphic import create_operator_html as create_html
import sys
import pickle

''' functions:

acc2homolog_list(accessionID, perc_ident, blastCacheFile, homologListFile)
    -> homologListFile with "accession" and "identity" keys

appendIntergenic(homologListFile) 
    -> homologListFile with "intergenic" and "regType" keys appended

appendOperatorMetadata(homologListFile, knownOperator)
    -> homologListFile with "operator", "invRepeat", and "deltaG" keys appended
'''

########  USER INPUTS  ########
acur_accession = "WP_011336736.1"
camr_accession = "BAA03510.1"
lrpr_accession = "WP_019744253"
alkx_accession = "AEM66515.1"
mmsr_accession = "AGI24556.1"
glpr_accession = "WP_157707983.1"
glprTHAF27 = "WP_152491583.1"
bm3r1_accession = "WP_013083972.1"

acc = bm3r1_accession
operatorFile = "bm3r1_all.txt"

    #having period in accession name screws things up
if acc[-2] == ".":
    acc = acc[0:-2]
print(acc)

perc_ident = 50

#Also, include a known operator for a regulator within this cluster within the /knownOperators directory

    #Use a known operator if there is one. Otherwise, just use the intergenic region

def operator_or_intergenic(acc,operatorFile):
    try:
        with open(f"knownOperators/{operatorFile}") as f:
            operator = f.read().replace('\n','')
            print("Known operator found: "+operator)
            return operator
    except:
        with open(f"cache/homolog_metadata/{acc}.pkl", mode="rb") as f:
            homologList = pickle.load(f)
            operator = homologList[0]["intergenic"]
            print(operator)
            print("No known operator found. Using intergenic region")
            return operator

def yes_or_no(question):
    reply = str(input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return True
    if reply[0] == 'n':
        return False
    else:
        return yes_or_no("Uhhh ... pls enter ")



appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
operator = operator_or_intergenic(acc,operatorFile)
appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator)
consensus_data = appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator)
create_html(consensus_data,"display/html_pages/"+acc)

"""
try:
    operator = operator_or_intergenic(acc,operatorFile)
    consensus_data = appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator)
    create_html(consensus_data,"display/html_pages/"+acc)
except:
    print("no intergenic data cached. Pulling intergenic data")
    try:
        yes_or_no('get intergenic region?')
        appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
        operator = operator_or_intergenic(acc,operatorFile)
        consensus_data = appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator)
        create_html(consensus_data,"display/html_pages/"+acc)
    except:
        print("no homolog data cached. Blasting protein")
        acc2homolog_list(acc, perc_ident)
        yes_or_no('get intergenic region?')
        appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
        operator = operator_or_intergenic(acc,operatorFile)
        consensus_data = appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator)
        create_html(consensus_data,"display/html_pages/"+acc)
"""



#######  PROGRAMS TO RUN ########
    #acc2homolog_list(accession, blastCacheFile, perc_ident, homologListFile)

    #appendIntergenic(homologListFile)
#appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
    
    #appendOperatorMetadata(homologListFile, knownOperator)
#appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", knownOperator)

