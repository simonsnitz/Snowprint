from acc2homologs import acc2homolog_list
from getIntergenic import appendIntergenic
from operator_utils import appendOperatorMetadata
import sys

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

acc = alkx_accession

    #having period in accession name screws things up
if acc[-2] == ".":
    acc = acc[0:-2]
print(acc)

perc_ident = 50

    #regulator_name = "lplr"

#Also, include a known operator for a regulator within this cluster within the /knownOperators directory



#######  FILES TO CREATE OR FETCH  ########
    #blastCacheFile = "cache/blastCache/"+regulator_name+str(perc_ident)+".xml"

    #homologListFile = "homolog_metadata/"+regulator_name+str(perc_ident)+".pkl"

#operator = regulator_name+".txt"
operator = "alkx.txt"

with open(f"knownOperators/{operator}") as f:
    knownOperator = f.read().replace('\n','')

print("Input operator: "+knownOperator)

#appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
#appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", knownOperator)


#"""
try:
    appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", knownOperator)
except:
    print("no intergenic data cached. Pulling intergenic data")
    try:
        appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
        appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", knownOperator)
    except:
        print("no homolog data cached. Blasting protein")
        acc2homolog_list(acc, perc_ident)
        appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
        appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", knownOperator)
#"""



#######  PROGRAMS TO RUN ########
    #acc2homolog_list(accession, blastCacheFile, perc_ident, homologListFile)

    #appendIntergenic(homologListFile)
#appendIntergenic(f"cache/homolog_metadata/{acc}.pkl")
    
    #appendOperatorMetadata(homologListFile, knownOperator)
#appendOperatorMetadata(f"cache/homolog_metadata/{acc}.pkl", knownOperator)

