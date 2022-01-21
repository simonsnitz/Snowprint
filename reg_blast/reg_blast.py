from processing.acc2homologs import acc2homolog_list
from processing.getIntergenic_batch import appendIntergenic_batch
from processing.operator_utils import appendOperatorMetadata

from display.append_operon import create_operon_html
from display.append_operator import create_operator_html

import imgkit
import webbrowser
from pathlib import Path
import pickle

    # Input parameters
tetr = "WP_000113282.1"
edcr = "AZI36850"
dszgr = "WP_058249973.1"
pnbx = "QWF22305.1"
sace_0303 = "CAL99652.1"
rosr = "WP_003861276.1"
testo = "EHN63435"
acc = testo

"""
    Input parameters:
"""
# regulators = [testo, tetr, sace_0303, rosr, dszgr, edcr, pnbx]
regulators = ["EHN63435", "CAM88408", "SACOL2593", "CAL99652.1", "WP_009944749.1", "WP_011014502.1", "WP_005379994.1", "WP_003416573.1", \
    "WP_011028825.1", "WP_010986683.1", "WP_003861276.1", "QWF22305.1", "WP_058249973.1", "WP_000802707.1", "AZI36850", "WP_012480398.1",\
    "WP_001295289.1", "NP_639852.1", "WP_087743244.1", "WP_011089760.1", "WP_003112936.1", "AAY80459", "AEM66515.1", "WP_009949582.1", \
    "WP_013366341.1", "WP_001807342.1", "WP_000113282.1", "WP_014859138", "WP_000113609.1", "WP_158491109", "WP_011336736.1", "WP_024719023.1",\
    "BAA03510.1", "NP_414847.3", "WP_169433046.1", "WP_010905492.1", "WP_011728885.1", "WP_011030045.1", "NP_252264.1", "NP_415533.1", \
    "WP_109753234.1", "WP_011015509.1", "AJS09379.1", "WP_011087095.1", "WP_052413891.1", "WP_011014162.1", "WP_004926705.1", "WP_011156222.1", \
    "BAK67179", "ACY33523", "WP_000805902.1", "WP_000224470.1", "WP_000762401.1", "WP_013636164.1", "WP_003514478.1", "WP_011291385.1", \
    "WP_011029024.1", "NP_388271.1", "WP_010920250.1", "WP_003243070.1", "WP_150675592.1", "WP_000174305.1", "WP_087081127.1", "CAG68545", \
    "WP_011731512.1", "WP_011594778.1", "ABP48117.1", "WP_087081125.1"]

regulators = [testo]

# Had an error caching blast results for WP_003416573.1, WP_087743244.1
# urllib.error.URLError: <urlopen error [Errno -3] Temporary failure in name resolution>

max_seqs = 50
perc_ident = 50
to_align = "find_inverted_repeat"
# Options: [ "find_inverted_repeat", "whole_interoperon", or custom string]


"""
ERRORS TO ADDRESS:

1. Had an error caching blast results for WP_003416573.1, WP_087743244.1
   urllib.error.URLError: <urlopen error [Errno -3] Temporary failure in name resolution>

2. 
"""


    # Having period in accession name screws things up. Remove it
if acc[-2] == ".":
    acc = acc[0:-2]

    # Cache file locations:
# BLAST: 'cache/blast_cache/{acc}.xml'
# Metadata: 'cache/homolog_metadata/{acc}.pkl'
# Updated Metadata: 'cache//homolog_metadata/updated_metadata/{acc}.pkl'
# Operator data: 'cache/operator_data/{acc}.pkl'
# HTML pages: 'cache/HTML/{date}.html


    # Takes a regulator accession ID and outputs a dictionary with predicted operator data
def acc2operator_data(acc, max_seqs):

    acc2homolog_list(acc, max_seqs)
    appendIntergenic_batch(acc)
    operator_data = appendOperatorMetadata(acc, to_align, perc_ident)
    return operator_data

    # Create operator object for each input regulator
operators = [acc2operator_data(reg, max_seqs) for reg in regulators]


    # Get operon object from cache
with open(f'cache/homolog_metadata/updated_metadata/{acc}.pkl', mode="rb") as f:
    operons = [i for i in pickle.load(f)]
operon = operons[0]
operon["organism"] = "placeholder"
operon = [operon]


    # update operon and operator JS files
create_operon_html(operon)
create_operator_html(operators) 


    # Open updated HTML page
p = Path("./display/Webpage")
fp = p / 'report.html'
html_name = str(fp.resolve())

webbrowser.open(html_name)


    # Create an image from HTML display
options = {
    'format': 'png',
    'width': 1000,
    'height':700,
    'encoding': "UTF-8"
}
imgkit.from_file(html_name, 'display/report.png', options=options)




#def create_display(operons, acc, consensus_data):
    #create_operon_html(operons,"display/html_pages/"+acc)
    #create_operator_html(consensus_data,"display/html_pages/"+acc+".html")
    #webbrowser.open('display/html_pages/'+acc+'.html')


### ALL INPUT PARAMETERS ###

    #BLAST
#hitlist_size = 100 (int input: default-100)
#perc_ident = [70,60] (int input: default-70)

    #ACCESSION
#acc = "WP_020114152.1" (string input)

    #OPERATOR
#option_1: provide predicted operator sequence (string input)
#option_2: look for inverted repeats (check box)
#option_3: align entire inter-operon region (check box)

    #DISPLAY
#show metadata (check box: default-yes)
#show operons (check box: default-yes)

    #EMAIL
#email address to send to (string input) --> form validation





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

RslR4 = "WP_020114152.1"
#acc = "WP_011014162.1"
