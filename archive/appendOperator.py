"""
***ABOUT***

    Prototype version of function used to extract an operator from an intergenic region by aligning a known operator of a homolog.
    
    An improved version of this function is found in the 'operator_utils.py' file

"""

import pickle
from  pprint import pprint
from Bio.pairwise2 import format_alignment, align


def findOperatorInUpstream(upstream, operator):
    """
        Do a local alignment with gap penalties:
            Identical characters add 2 points
            Non-Identical characters deduct 1 point
            Opening Gap deduct 1 point
            Extending a Gap  deduct 1 point
        Arguments:
            upstream: SeqRecord of upstream DNA region
            operator: SeqRecord of operator DNA sequence
        Return:
            Void
    """
    #Decides if input is string or SeqRecord
    if type(operator) == str:
        pass
    else:
        operator = operator.seq

    operator_length = len(operator)

    #this function is needed to extract the ENTIRE aligned region.
        #Otherwise, mismatched ends will not be extracted.
    def extractOperator(upstream, operator):
        beginning = 0
        for i in operator:
            if i == '-':
                beginning += 1
            else:
                break
        end = (beginning + operator_length)
        return upstream[beginning:end]


    try:
        upstr_align, op_align, score, startPos, endPos = \
        align.localms(upstream, operator, 2, -0.5, -100, 0)[0]
    except:
        return "no intergenic region found"

    #Heavily penalizing gap opening avoids errors with downstream data processing,
    #but may miss out on interesting biological features

    #returns the aligned operator sequence if a similarity threshold is met.
        #score threshold should be tuned.

    if score > 7:
        operator = extractOperator(upstr_align, op_align)
        return operator, score
    else:
        #print('too low', score)
        return "Score is too low: "+str(score)


if __name__ == "__main__":

    with open('homolog_list2.pkl', mode="rb") as f:
        homologList = pickle.load(f)

    operator = "GTATATCGCAGATATAG"

    for i in homologList:
        i["operator"] = findOperatorInUpstream(i["intergenic"], operator)
        print(i["operator"])
    
    maxScore = len(operator)*2
    print("max score = "+str(maxScore))
    #pprint(homologList)
