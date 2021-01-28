from processing.acc2homologs import acc2homolog_list
from processing.getIntergenic import appendIntergenic
from processing.operator_utils import appendOperatorMetadata

from display.render_operator_graphic import create_operator_html
from display.render_operon_graphic import create_operon_html

import sys
import pickle
import platform


####### DEFINE BLAST/FILTER PARAMETERS #########
hitlist_size = 100
perc_ident = 70


########  USER INPUTS  ########
ramr = "WP_000113609.1"
acur = "WP_011336736.1"
camr = "BAA03510.1"
lplr = "WP_019744253"
alkx = "AEM66515.1"
bm3r1 = "WP_013083972.1"
sco4850 = "WP_011029905.1"
mmsr = "AGI24556.1"
hpdr = "WP_010951464.1"
alku = "WP_009475225.1"
pocr = "WP_000622326.1"
glpr = "WP_157707983.1"
glprTHAF27 = "WP_152491583.1"
glprPAED = "WP_165195465.1"
glprTHAS = "WP_106472839.1"

#ProXes
ttgr = "WP_010952495.1"
cgmr = "WP_011015249.1"
acrr = "WP_000101737.1"
mtrr = "WP_003693763.1"
eilr = "WP_013366341.1"
#cgmr = "WP_011015249.1"
ebrr = "WP_003976902.1"
smet = "WP_005414519.1"

acc = smet


###### OPERATOR INPUT #######
operatorFile = "camr.txt"

with open(f"knownOperators/{operatorFile}") as f:
    operator = f.read().replace('\n', '')

#having period in accession name screws things up. Remove it
if acc[-2] == ".":
    acc = acc[0:-2]
print(acc)


##### FUNCTIONS #####


#acc2homolog_list(acc, hitlist_size)
#appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")

consensus_data = appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator, perc_ident) 

print(consensus_data)

#create_operator_html(consensus_data,"display/html_pages/"+acc+".html")


'''
if platform.release() == "4.4.0-18362-Microsoft":
    create_operon_html(operons,"../../../../mnt/c/Users/simon.doelsnitz/"+acc)
    create_operator_html(consensus_data,"../../../../mnt/c/Users/simon.doelsnitz/"+acc+".html")
else:
    create_operon_html(operons,"display/html_pages/"+acc)
    create_operator_html(consensus_data,"display/html_pages/"+acc+".html")
'''




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
