"""
***DEPRECATED***

SEARCH_TO_GRAPHIC replaced this constructor script

"""

#TODO:

#INPUT: a chemical (SMILES or InChIKey format)

#OUTPUT: a list of candidate sensors for that chemical, and their corresponding regulators

'''
chemical2enzymes.py         -get list of enzymes that act on input chemical
->
acc2operon.py               -find regulators associated with each enzyme (if any)
->
getIntergenic.py            -pull out intergenic region for each regulator
->
#####################################################################################   optional. For better operator prediction
acc2homologs.py             -get regulator homologs (accession ID and % identity)
->
getIntergenic.py            -get intergenic region for each homolog
->
#####################################################################################
operator_utils.py           -get predicted operator for enzyme-associated regulator

===> Return list of enzyme:regulator:operator sets (with % confidence in operator prediction)
