from acc2homologs import acc2homolog_list
from getIntergenic import appendIntergenic
from operator_utils import appendOperatorMetadata

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

accession = acur_accession

perc_ident = 50

regulator_name = "camr"

#Also, include a known operator for a regulator within this cluster within the /knownOperators directory



#######  FILES TO CREATE OR FETCH  ########
#blastCacheFile = "cache/blastCache/"+regulator_name+str(perc_ident)+".xml"

homologListFile = "homolog_metadata/"+regulator_name+str(perc_ident)+".pkl"

operator = regulator_name+".txt"
with open(f"knownOperators/{operator}") as f:
    knownOperator = f.read().replace('\n','')



#######  PROGRAMS TO RUN ########
#acc2homolog_list(accession, blastCacheFile, perc_ident, homologListFile)

appendIntergenic(homologListFile)

#appendOperatorMetadata(homologListFile, knownOperator)

