import sys
import os
import os.path
import shutil

import jellyfish


def compare(codon1, codon2):
    s = jellyfish.hamming_distance(codon1, codon2)
    return s


synCodons = {
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'Q': ['CAA', 'CAG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCT', 'CCG', 'CCA', 'CCC'],
    'K': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'T': ['ACC', 'ACA', 'ACG', 'ACT'],
    'F': ['TTT', 'TTC'],
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'G': ['GGT', 'GGG', 'GGA', 'GGC'],
    'I': ['ATC', 'ATA', 'ATT'],
    'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'H': ['CAT', 'CAC'],
    'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'W': ['TGG'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'E': ['GAG', 'GAA'],
    'Y': ['TAT', 'TAC']}


d = os.getcwd()
os.chdir("..")


#works on Linux, but not windows. Uggh!
sys.path.insert(1, '/home/simon/git/tuSeek')

from acc2homologs import homologs2residueFrequency as a2r

#ramr
acc = "WP_000113609"

#camr
#acc = "BAA03510"

#seq = a2h.accID2sequence(acc)
#print(seq)

residue_freq = a2r(acc,70)


res = list(residue_freq[30].keys())


def residues2libCodon(residues):

    #get synonomous codons
    codons = []
    for i in residues:
        codons.append(synCodons[i])

    #initiate codon comparison (compares one codon to one other codon)
    codonPairs = []
    for i in codons[0]:
        for j in codons[1]:
            if compare(i,j) == 1:
                codonPairs.append([i,j])
            #Need to catch exception where lowest hamming distance is more than 1

    #continue codon comparison (compares one codon to an array of codons)  Make this a function.
    allPairs = []
    for i in codonPairs:
        for k in codons[2]:
            score = 0
            for j in i:
                score += compare(j,k)
            allPairs.append([score,k,i])

    lowestScore = min([ x[0] for x in allPairs])
        #filter codon pairs with the lowest scores.
    bestPairs = [ x for x in allPairs if x[0] == lowestScore ]
    #Cant append within a list comprehension?


    print(bestPairs)

residues2libCodon(res)

'''
minScore = min(scores)

bestScores = []
for i in codons[0]:
    for j in codons[1]:
        if compare(i,j) == minScore:
            bestScores.append([i,j])

#print(bestScores)
'''


