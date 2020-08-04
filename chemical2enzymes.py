#!/usr/bin/env python
import Bio
import time
from Bio.KEGG.REST import *
from Bio.KEGG import REST
from Bio.KEGG import Enzyme
import re
from Bio import Entrez, SeqIO

Entrez.email = "brianna.flynn@utexas.edu"

# Experimenting with using KEGG database to pull out EC numbers
'''
main function: Chemical2Enzymes(input_chemical: str, filter criteria: str)

available functions:                                                                                                                                                                                    

chemicalInput(input_chemical : str) 
EC_converter(EC_numbers (list)) # formats EC numbers to query NCBI
ID_dictionaries(EC_IDs) # takes formatted EC nums, returns gene IDs
Gene_Dict(IDs) # takes gene IDs, returns a dictionary with their accession numbers as keys and taxonomies as values
Filter(search_key: str, gene_dict) # takes a taxonomic criteria (bacteria, fungi, etc) and a dict of {accession: taxonomy} and returns the list of accessions as well as a filtered list of accessions \
based on the search criteria

'''
def chemicalInput(input_chemical: str):
    enzymes = kegg_list("enzyme", org=None).read()
    
    enz_descrip = []
    enz_entry = []

    for line in enzymes.rstrip().split("\n"):
        entry, description = line.split("\t")
        if input_chemical in description:
            enz_descrip.append(description)
            enz_entry.append(entry)
    print("Input Chemical : " + str(input_chemical))        
    return enz_entry


# Looking for all enzyme entries in Kegg that relate to geraniol
# Index the  second entry from the list of enzyme ID's and append it to second_entry variable in order to convert EC number to NCBI friendly format
# need to fix this if there is just one number, not list

def EC_converter(EC_num_list):
    EC_str = '[EC]'
    EC_num = []
    EC_IDs = []
    for entry in EC_num_list:
        EC_num.append(entry[3:])
        
    for num in EC_num:
        EC_IDs.append(num+EC_str)
    
    return EC_IDs

# This takes the EC_IDs from kegg and converts them into a format that can be used to query the ncbi protein db for the EC number's corresponding genes

def ID_dictionaries(EC_IDs):
    ID_dicts=[]
    for id in EC_IDs:
        handle=Entrez.esearch(db='protein', term=id)
        entries=Entrez.read(handle)
        ID_dicts.append(entries)
        handle.close()

    IDs=[]
    for dict in ID_dicts:
        for key, value in dict.items():
            if key == 'IdList':
                IDs.append(value)
    IDs = [y for x in IDs for y in x]
    
    return IDs

# Takes the EC_IDs of the chemical, returns the IDs of genes associated with the EC_IDs obtained from kegg

def GeneDict(IDs):    
    Genes=[]
    accessions=[]
    taxonomies=[]
    
    for id in IDs:
        handle = Entrez.efetch(db='protein', id=id, retmode='xml')
        records = Entrez.read(handle)
        Genes.append(records)
        handle.close()    
    Genes=[y for x in Genes for y in x]
    
    for gene in Genes:
        for k, v in gene.items():
            if k == 'GBSeq_primary-accession':
                accessions.append(v)
    for gene in Genes:
        for k, v in gene.items():
            if k == 'GBSeq_taxonomy':
                taxonomies.append(v)
    gene_dict = {k: v for k, v in zip(accessions, taxonomies)}
    
    return gene_dict

# GeneDict returns a dictionary with accession IDs of each gene pulled from ID dictionary as key, and the taxonomy of the gene as the value. This is used \
# as input to filter accession IDs based on what org enzyme comes from

def Filter(search_key : str, gene_dict):
    
    query = search_key
    all_acc = [key for key, val in gene_dict.items()]
    filtered = [key for key, val in gene_dict.items() if query in val]
    print("")
    print("All Enzyme Accession IDs found: " + str(all_acc))
    print("")
    print("Filter criteria : " + str(query))
    print("Filtered Accession IDs : "  + str(filtered))

    return filtered


def Chemical2Enzymes(input_chemical : str, filter_criteria: str):
    Chemical = chemicalInput(input_chemical)
    ECs = EC_converter(Chemical)
    Gene_Ids = ID_dictionaries(ECs)
    Gene_Dicts = GeneDict(Gene_Ids)
    Filter(filter_criteria, Gene_Dicts)

# This filters the GeneDict based on the criteria specified as the search key. Prints both the accession IDs of the non-filtered and filtered GeneDict, and then returns the filtered list 

if __name__ == "__main__":
    
    Test = Chemical2Enzymes("geraniol", "Bacteria")
    print(Test)
    
# TO DO: Test this with mutliple chemicals. 
# Make edits to account for different scenarios (only one EC_ID, program fails to return anything for specified chemical, no bacterial enzymes are present, etc)

