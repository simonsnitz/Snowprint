from acc2homologs import acc2homolog_list
from getIntergenic import appendIntergenic
from operator_utils import appendOperatorMetadata

from display.templates.render_operator_graphic import create_operator_html
from display.templates.render_operon_graphic import create_operon_html
import sys
import pickle
import platform


####### DEFINE BLAST/FILTER PARAMETERS #########
hitlist_size = 100
perc_ident = [90,70,60]


########  USER INPUTS  ########
acur = "WP_011336736.1"
camr = "BAA03510.1"
lplr = "WP_019744253"
alkx = "AEM66515.1"
bm3r1 = "WP_013083972.1"
sco4850 = "WP_011029905.1"
mmsr = "AGI24556.1"
alku = "WP_009475225.1"
pocr = "WP_000622326.1"
glpr = "WP_157707983.1"
glprTHAF27 = "WP_152491583.1"
glprPAED = "WP_165195465.1"
glprTHAS = "WP_106472839.1"

acc = mmsr


###### OPERATOR INPUT #######
operatorFile = "mmsr.txt"
#operatorFile = "None"


#having period in accession name screws things up. Remove it
if acc[-2] == ".":
    acc = acc[0:-2]
print(acc)


##### FUNCTIONS #####
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



#acc2homolog_list(acc, hitlist_size)
appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")

operator = operator_or_none(acc,operatorFile)

consensus_data = [ appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator, i) 
    for i in perc_ident
    ]

with open(f"cache/homolog_metadata/{acc}.pkl", mode='rb') as f:
    data = pickle.load(f)
    operons = data

if platform.release() == "4.4.0-18362-Microsoft":
    create_operon_html(operons,"../../../../mnt/c/Users/simon.doelsnitz/"+acc)
    create_operator_html(consensus_data,"../../../../mnt/c/Users/simon.doelsnitz/"+acc+".html")
else:
    create_operon_html(operons,"display/html_pages/"+acc)
    create_operator_html(consensus_data,"display/html_pages/"+acc+".html")





"""   Automate workflow. Check if things are cached. If not, get the necessary data. Fails if there are any bugs.
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
def yes_or_no(question):
    reply = str(input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return True
    if reply[0] == 'n':
        return False
    else:
        return yes_or_no("Uhhh ... pls enter ")

