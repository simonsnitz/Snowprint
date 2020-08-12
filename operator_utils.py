import pickle
from  pprint import pprint
from Bio.pairwise2 import format_alignment, align
import primer3
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from io import StringIO




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
        #Set score cutoff to be 40% of max. Arbitrary, but seems reasonable.
    max_score = 2*operator_length
    score_cutoff = max_score*0.4

    if score > score_cutoff:
        operator = extractOperator(upstr_align, op_align)
        return [operator, score]
    else:
        return [None,0]


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


def getBestInvertedRepeat(operator):
    stringency = [10,9,8,7,6,5,4]
    if operator != None:
        for i in stringency:
            invRepeat = findInvertedRepeat(operator, i)
            if invRepeat:
                return [invRepeat, i]

'''
def alignIntergenic(homologList):
    allIntergenic = []
    for i in homologList:
        if "intergenic" in i.keys():
            if i["identity"] >= 50:
                #print(i["intergenic"], i["identity"])
                allIntergenic.append(SeqRecord(Seq(i["intergenic"], generic_dna), id=str(i["identity"])))

    SeqIO.write(allIntergenic, "cache/camrIntergenic.fasta", "fasta")

    cline = MuscleCommandline(input="cache/camrIntergenic.fasta", gapopen=-200.0, gapextend=-10.0)
    std_out, std_err = cline()
    align = AlignIO.read(StringIO(std_out), "fasta")
    print(align[2].seq)
'''
            #filter by identity or by alignment score? Both? Some alignments from ~80% homologs have crap scores
def getConsensus(homologList, max_ident=100, min_ident=70):      #DUUUDE! Mad list comprehension!!!
    
    allOperators = [ i["operator"] for i in homologList
            if "operator" in i.keys() and i["identity"] >= min_ident and i["identity"] <= max_ident and i["score"] != 0]

    #bases_at_each_position = [ [i[pos] for i in allOperators if len(i) == 27]   #added if statement do deal with weirdness
    bases_at_each_position = [ [i[pos] for i in allOperators]   #added if statement do deal with weirdness
            for pos in range(0, len(allOperators[0])) ] 

    def mostFrequent(List):         #function to make it more readable
        return max(set(List), key=List.count)

    consensusOperator = "".join(mostFrequent(i) for i in bases_at_each_position)
    print(bases_at_each_position)

    return consensusOperator
    #return getBestInvertedRepeat(consensusOperator)


def getMostSymmetric(homologList, min_identity=65):

    bestIR = max([i["invRepeat"][1] for i in homologList if i["invRepeat"] != None if "invRepeat" in i.keys()  if i["identity"] >= min_identity])

    for i in homologList:
        if i["invRepeat"] != None:
            if i["invRepeat"][1] == bestIR:
                mostSymmetric = i["invRepeat"]
                break

    return mostSymmetric


def appendOperatorMetadata(homologListFile, knownOperator):
    
    with open(f'{homologListFile}', mode="rb") as f:
        homologList = pickle.load(f)
        for i in homologList:
            #try:
            i["operator"], i["score"] = findOperatorInIntergenic(i["intergenic"], knownOperator)
            i["invRepeat"] = getBestInvertedRepeat(i["operator"])
            
            #if i["score"] != 0:
                #print(i["operator"])
                #print(i["identity"])
                #print(i["invRepeat"])
       
        
        print("      consensus sequence")
        print(getConsensus(homologList))
        return homologList
            

if __name__ == "__main__":

    ramr_operator = "ATGAGTGAGTACGCAACTCAT"
    camr_operator = "GTATATCGCAGATATAG"

    #print('getting operators')

    #homologList = appendOperatorMetadata('homolog_metadata/ramr50.pkl', ramr_operator)
    with open('homolog_metadata/ramr50.pkl', mode='rb') as f:
        homologList = pickle.load(f)
    
    #with open('ramr50.pkl', mode="wb") as f:
    #    pickle.dump(homologList, f)
   
    #pprint(homologList)
    
    consensus = getConsensus(homologList)
    print("consensus operator: "+str(consensus))

    mostSymmetric = getMostSymmetric(homologList)
    print("most symmetric operator: "+str(mostSymmetric))
    
    #alignIntergenic(homologList)
    
    #print("max score = "+str(len(operator)*2)
    
