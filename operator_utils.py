import pickle
from  pprint import pprint
from Bio.pairwise2 import format_alignment, align
import primer3

def findOperatorInIntergenic(intergenic, operator):
    """
        Do a local alignment with gap penalties:
            Identical characters add 2 points
            Non-Identical characters deduct 0.5 point
            Opening Gap deduct 100 points
            Extending a Gap  deduct 0 points
        Arguments:
            intergenic: Intergenic DNA region (string)
            operator: Known operator for regulator or its homolog (string)
    """
    operator_length = len(operator)
    '''
        This function is needed to extract the ENTIRE aligned region.
        Otherwise, mismatched ends will not be extracted.
    '''
    def extractOperator(intergenic, operator):
        beginning = 0
        for i in operator:
            if i == '-':
                beginning += 1
            else:
                break
        end = (beginning + operator_length)
            #can change position of output sequence to include more or less
        return intergenic[beginning-3:end+3]

    try:
        upstr_align, op_align, score, startPos, endPos = \
        align.localms(intergenic, operator, 2, -0.5, -100, 0)[0]
    except:
        return "no intergenic region found"
    '''
    Heavily penalizing gap opening avoids errors with downstream data processing, but may miss out on interesting biological features
    Returns the aligned operator sequence if a similarity threshold is met. Score threshold (7) should be tuned.
    '''
    if score > 7:
        operator = extractOperator(upstr_align, op_align)
        return operator, score
    else:
        #print('too low', score)
        return "Score is too low: "+str(score)


def complement(sequence):
    compDNA = {"A":"T",
            "C":"G",
            "T":"A",
            "G":"C"}
    complement = ""
    for i in sequence:
        try:
            complement += compDNA[i]
        except:
            print('non standard base found!!!')
            break
    return complement


def findInvertedRepeat(intergenic, size):
    for i in range(0,len(intergenic)-((2*size))):
        seq = intergenic[i:i+size]
        revCompSeq = complement(seq)[::-1]
        maxAllowedSpacer = i+size+8
            
        for j in range(i+size,maxAllowedSpacer):
            compare = intergenic[j:j+size]
            if compare == revCompSeq:
                
                beginning = intergenic[0:i].lower()
                seqF = seq.upper()
                middle = intergenic[i+size:j].lower()
                seqR = compare.upper()
                end = intergenic[j+size:].lower()
                operator = beginning+seqF+middle+seqR+end
                return operator


def getBestInvertedRepeat(intergenic):
    stringency = [10,9,8,7,6,5,4]
    if intergenic != None:
        for i in stringency:
            operator = findInvertedRepeat(intergenic, i)
            if operator:
                return [operator, i]


def appendOperatorMetadata(homologListFile, knownOperator):
    
    with open(f'{homologListFile}', mode="rb") as f:
        homologList = pickle.load(f)
        for i in homologList:
            try:
                i["operator"] = findOperatorInIntergenic(i["intergenic"], knownOperator)
                i["invRepeat"] = getBestInvertedRepeat(i["operator"][0])
                    #Calculate the free energy (deltaG) of hairpin formation within the operator sequence.
                i["deltaG"] = primer3.calcHairpin(i["invRepeat"][0]).dg
                print(i)
            except:
                print("data not found")
        
        return homologList
            

if __name__ == "__main__":

    with open('homolog_ramrIntergenic.pkl', mode="rb") as f:
        homologList = pickle.load(f)

    #camr
    #operator = "GTATATCGCAGATATAG"

    #ramr
    operator = "ATGAGTGAGTACGCAACTCAT"

    for i in homologList:
        try:
            i["operator"] = findOperatorInUpstream(i["intergenic"], operator)
            i["invRepeat"] = getBestOperator(i["operator"][0])
            i["deltaG"] = calcHairpin(i["invRepeat"][0])
            #pprint(i)
            print(i["identity"],i["regType"], i["invRepeat"], i["deltaG"])
        except:
            print("data not found")
    
    maxScore = len(operator)*2
    print("max score = "+str(maxScore))
    
    with open('ramr_homologOperators.pkl', mode="wb") as f:
        pickle.dump(homologList, f)
