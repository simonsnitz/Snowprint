
#TODO:

#INPUT: (1) an enzyme accession ID      (2) % identity cutoff to search for enzyme homologs
#Rationale: if no regulator associated with target enzyme, maybe there is one associated with a closely related homolog

#OUTPUT: a list of candidate sensors associated with that enzyme / its homologs, and their corresponding regulators

'''
acc2homologs.py             -get list of enzyme homologs (supposedly perform same chemistry)
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
'''

from acc2homologs import acc2homolog_list

from acc2operon import enzyme_acc2regulators




