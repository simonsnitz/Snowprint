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
perc_ident = [70,60]


########  USER INPUTS  ########
tetr = "WP_000113282.1"
acur = "WP_011336736.1"
camr = "BAA03510.1"
lplr = "WP_019744253"
alkx = "AEM66515.1"
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
ramr = "WP_000113609.1"
qacr = "WP_001807342.1"
eilr = "WP_013366341.1"
bm3r1 = "WP_013083972.1"
ttgr = "WP_010952495.1"
lmrr = "WP_011834386.1"
ebrr = "WP_003976902.1"
mexr = "WP_003114897.1"
#ladr = "WP_003721913.1"
#vcer = "WP_001264144.1"
mtrr = "WP_003693763.1"
acrr = "WP_000101737.1"
mepr = "WP_000397416.1" #nonetype not iterable for regulator 95...
smet = "WP_005414519.1" #homologs and intergenic cached
nald = "WP_003092152.1"
lfrr = "WP_003897643.1" #homologs and intergenic cached, error with operator
mdtr = "WP_003242592.1"
mexr = "WP_003088626.1"
aden = "AGV28567.1"
nfxb = "WP_033999467.1" #data collected. No operator found.
sco4008 = "WP_011029378.1"
rv3066 = "WP_003416005.1" #interested pattern found. compare operator to reference.
vexr = "WP_001884097.1" #homologs and intergenic cached. No operator found
cgmr = "WP_011015249.1" #homologs and intergenic cached. No operator found
rv0302 = "WP_003401571.1" #homologs and intergenic cached. Error with outputting HTML file. one intergenic region REALLY long
bdtr = "WP_011089760.1" #looks like a pretty decent operator prediction
rv1353c = "WP_003898836.1" #nothing convincing

oqxr = "WP_000888203.1" #homologs and intergenic cached. Error with outputting HTML file. one intergenic region REALLY long
smvr = "WP_004888643.1" #small IR found. Map to known operator
bper = "WP_004526226.1" #perfect operator found!

mgra = "WP_001283444.1" #intergenic cached. No predicted operators look convinving
mexl = "WP_003092468.1" #perfect operator found
satr = "WP_012027921.1" #perfect operator found. Agrees with literature :)
cmer = "WP_002857627.1" #perfect opertator found
fepr = "AHN05314.1" #perfect operator found!
ttgt = "WP_012052586.1" #intergenic cached. No consensus operator found.
ttgv = "WP_014003968.1" #a short 4-bp operator found. Follow up
rv0678 = "WP_003403442.1" #decent 4bp IR found
brta = "WP_003723289.1"
lmo1618 = "WP_003723583.1"
brer = "WP_000777169.1" #data collected. operator predictor didn't pick up anything

lmra = "WP_003246449.1" #no great operator found. homologs/intergenic region cached.
acrs = "WP_001129518.1" #IR found but it may be a terminator. map to known sequence.
bce2991 = "WP_001107304.1" #perfect IR found. Same as Bm3R1.
abacrr = "WP_004929123.1"
tm1030 = "WP_010865247.1"
stacrr = "WP_000101755.1"
cac3606 = "WP_010966869.1" #great IR found!
acrr1 = "WP_010880658.1" #not enough close homologs. only 3 seqs aligned
ybih = "WP_001025254.1" #long IR found, but it's not conserved.
smu134 = "NP_720607.1" #small IR found, decent but not excellent conservation
sco0520 = "WP_003978343.1"
tcar = "WP_001832914.1"
sco4122 = "WP_003974850.1"

acc = sco4122


###### OPERATOR INPUT #######
#operatorFile = "hpdr3.txt"
operatorFile = "None"


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
#appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")

#operator = operator_or_none(acc,operatorFile)

operator = "CACCTTCGAACTTTAGCTTCTAAGTCTT"

consensus_data = [ appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", operator, i) 
    for i in perc_ident
    ]

with open(f"cache/homolog_metadata/{acc}.pkl", mode='rb') as f:
    data = pickle.load(f)
    operons = data

create_operon_html(operons,"display/html_pages/"+acc)
create_operator_html(consensus_data,"display/html_pages/"+acc+".html")


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
