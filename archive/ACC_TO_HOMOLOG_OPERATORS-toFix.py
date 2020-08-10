"""
***DEPRECATED***

    OPTIMIZE_OPERATOR function already does this. It's redundant to include this as a separate file

"""


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

#accession = "BAA03510.1"
#acur
accession = "WP_011336736.1"

perc_ident = 50

    #file must be xml format
blastCacheFile = 'cache/blastCache/acur50.xml'

    #file must be pkl format
homologListFile = 'homolog_metadata/acur.pkl'

operator = "acur.txt"
with open(f"knownOperators/{operator}") as f:
    knownOperator = f.read().replace('\n','')


#acc2homolog_list(accession, blastCacheFile, perc_ident, homologListFile)
acc2homolog_list(accession, perc_ident)

appendIntergenic(homologListFile)

appendOperatorMetadata(homologListFile, knownOperator)
