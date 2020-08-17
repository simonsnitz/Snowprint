from acc2homologs import acc2homolog_list
from getIntergenic import appendIntergenic
from operator_utils import appendOperatorMetadata
from display.operator_graphic import create_operator_html as create_html
import sys
import pickle
import platform

''' functions:

acc2homolog_list(accessionID, perc_ident, blastCacheFile, homologListFile)
    -> homologListFile with "accession" and "identity" keys

appendIntergenic(homologListFile) 
    -> homologListFile with "intergenic" and "regType" keys appended

appendOperatorMetadata(homologListFile, knownOperator)
    -> homologListFile with "operator", "invRepeat", and "deltaG" keys appended
'''

########  USER INPUTS  ########
acur = "WP_011336736.1"
camr = "BAA03510.1"
lrpr = "WP_019744253"
alkx = "AEM66515.1"
bm3r1 = "WP_013083972.1"
sco4850 = "WP_011029905.1"
mmsr = "AGI24556.1"

glpr = "WP_157707983.1"
glprTHAF27 = "WP_152491583.1"
glprPAED = "WP_165195465.1"
glprTHAS = "WP_106472839.1"

acc = sco4850

#operatorFile = "mmsr.txt"
operatorFile = "None"

    #having period in accession name screws things up
if acc[-2] == ".":
    acc = acc[0:-2]
print(acc)

hitsize_list = 100
perc_ident = 90

#Also, include a known operator for a regulator within this cluster within the /knownOperators directory

    #Use a known operator if there is one. Otherwise, just use the intergenic region

def operator_or_none(acc,operatorFile):
    try:
        with open(f"knownOperators/{operatorFile}") as f:
            operator = f.read().replace('\n','')
            print("Known operator found: "+operator)
            return operator
    except:
        operator = None
        print("no known operator. Finding inverted repeats")
        return operator

def yes_or_no(question):
    reply = str(input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return True
    if reply[0] == 'n':
        return False
    else:
        return yes_or_no("Uhhh ... pls enter ")



#acc2homolog_list(acc, hitsize_list)
#appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
operator = operator_or_none(acc,operatorFile)
consensus_data = appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator, perc_ident)
#print(consensus_data)

if platform.release() == "4.4.0-18362-Microsoft":
    create_html(consensus_data,"../../../../mnt/c/Users/simon.doelsnitz/"+acc)
else:
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

