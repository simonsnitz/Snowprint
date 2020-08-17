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

#TODO

#score_cutoff is arbitrarily set to 40% of max score. Is that relevant? (line #56)




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



#Functions for finding inverted repeats within a given sequence

def findInvertedRepeat(intergenic, size):
    operators = []
    for i in range(0,len(intergenic)-((2*size))):
        seq = intergenic[i:i+size]
        revCompSeq = complement(seq)[::-1]
        maxAllowedSpacer = i+size+8

        for j in range(i+size,maxAllowedSpacer):
            compare = intergenic[j:j+size]
            if compare == revCompSeq and i>6:   #added filter so operator isn't at very edge of interoperon region
                
                #'''   #for formatting a string, with bases in inverted repeat upper case
                beginning = intergenic[i-5:i].lower()
                seqF = seq.upper()
                middle = intergenic[i+size:j].lower()
                seqR = compare.upper()
                end = intergenic[j+size:j+size+5].lower()
                print(beginning+seqF+middle+seqR+end)
                #'''
                try:
                    beginning = intergenic[i-10:i].lower()
                except:
                    beginning = intergenic[i-5:i].lower()
                seqF = seq.upper()
                middle = intergenic[i+size:j].lower()
                seqR = compare.upper()
                try:
                    end = intergenic[j+size:j+size+10].lower()
                except:
                    end = intergenic[j+size:j+size+5].lower()
                operator = beginning+seqF+middle+seqR+end

                operators.append(operator)
    return operators


def getBestInvertedRepeat(operator):
    stringency = [10,9,8,7,6,5,4]
    if operator != None:
        for i in stringency:
            invRepeat = findInvertedRepeat(operator, i)
            if invRepeat:
                return [invRepeat, i]

#End of inverted repeat functions

'''     
	#Alignment function using MUSCLE. Not working so well for operator prediction...
def alignIntergenic(homologList):
    allIntergenic = []
    allIdentities = []
    for i in homologList:
        if "intergenic" in i.keys():
            if i["identity"] >= 50:
                #print(i["intergenic"], i["identity"])
                allIntergenic.append(SeqRecord(Seq(i["intergenic"], generic_dna), id=str(i["identity"])))
                allIdentities.append(i["identity"])

    SeqIO.write(allIntergenic, "cache/intergenic.fasta", "fasta")

    cline = MuscleCommandline(input="cache/intergenic.fasta", gapopen=-100000.0, gapextend=-1000.0)
    std_out, std_err = cline()
    align = AlignIO.read(StringIO(std_out), "fasta")
        #formatting into list of dictionaries such that it can be used by the "get consensus" function
    msaHomologs = [ {"operator": align[i].seq, "identity": allIdentities[i], "score":i+1}
            for i in range(0,len(align))
            ]
    for i in msaHomologs:
        print(i["operator"])
    return msaHomologs
'''


def getConsensus(homologList, max_ident=100, min_ident=70):      #DUUUDE! Mad list comprehension!!!
    
            #filter by identity and by alignment score Some alignments from ~80% homologs have crap scores
    allOperators = [ i["operator"] for i in homologList
            if "operator" in i.keys() and i["identity"] >= min_ident and i["identity"] <= max_ident and i["score"] != 0]

    print("allOperators: "+str(allOperators))

	#initialize list of dictionaries
    baep = [{base:1}
            for base in allOperators[0] 
        ]

	#populate dataset with base representation from all input operators
    for operator in allOperators[1:]:
        if len(operator) == len(allOperators[0]):
            for pos in range(0, len(operator)):
                base = operator[pos]
                #if base != '-':
                try:
                    baep[pos][base] +=1
                except:
                    baep[pos].update({base:1})


    max_values = [max(baep[pos].values()) for pos in range(0,len(baep))]
    max_score = max(max_values)
	#convert base conservation scores as a percent of max
    max_values_percent = [ round(i/max_score,2) for i in max_values ] 

    def get_key(my_dict,val):
        for key, value in my_dict.items():
            if val == value:
                return key

        return "key doesn't exist"

	#create list of most conserved bases at each position
    consensusSeq = [ get_key(baep[pos], max_values[pos])
        for pos in range(0, len(allOperators[0]))
    ]

	#dictionary containing the base and it's score at each position
    consensus_data = [{"base":consensusSeq[i] , "score":max_values_percent[i]} 
            for i in range(0,len(max_values))
        ]


	#separate code block to format consensus sequence. If base highly conserved, it gets capitalized
    def ifConserved(consensus,max_score):
        if consensus["score"] == max_score:
            return consensus["base"].upper()
        else:
            return consensus["base"].lower()

    print("max score is: "+str(max_score))
    
    formattedSeq = "".join( ifConserved(consensus_data[pos], max_score) for pos in range(0,len(max_values)))
    print(formattedSeq)
	#end of formatting consensus seq function


    return consensus_data



def getMostSymmetric(homologList, min_identity=50):

    bestIR = max([i["invRepeat"][1] for i in homologList if i["invRepeat"] != None if "invRepeat" in i.keys()  if i["identity"] >= min_identity])

    for i in homologList:
        if i["invRepeat"] != None:
            if i["invRepeat"][1] == bestIR:
                mostSymmetric = i["invRepeat"]
                break

    return mostSymmetric


def appendOperatorMetadata(homologListFile, knownOperator, perc_ident):
    
    with open(f'{homologListFile}', mode="rb") as f:
        homologList = pickle.load(f)
        
        if knownOperator == None:
            operators = getBestInvertedRepeat(homologList[0]["intergenic"])[0]
            print("inverted repeat found: "+ str(operators))
        else:
            operators = [knownOperator]
        
        consensus_data = []
        for operator in operators:
            for i in homologList:
                i["operator"], i["score"] =  findOperatorInIntergenic(i["intergenic"], operator.upper())
            consensus_data.append(getConsensus(homologList, min_ident=perc_ident))
			
        
        return consensus_data
        
        ''' #for single operator is known
        for i in homologList:
            i["operator"], i["score"] = findOperatorInIntergenic(i["intergenic"], operator)
            i["invRepeat"] = getBestInvertedRepeat(i["operator"])
        consensus_data = getConsensus(homologList)
        print("got operator, inverted repeat, and consensus data. Sending out consensus data")
        return consensus_data
        '''
    

if __name__ == "__main__":

    ramr_operator = "ATGAGTGAGTACGCAACTCAT"
    camr_operator = "GTATATCGCAGATATAG"
    

    acc_operator = "bm3r1.txt"
    #with open(f'knownOperators/{acc_operator}', mode='r') as f:
    #    operator = f.read()

    bm3r1 = "WP_013083972.pkl"
    camr = "BAA03510.pkl"
    acc = camr
    with open(f'cache/homolog_metadata/{acc}', mode='rb') as f:
        homologs = pickle.load(f)
        
        operators = getBestInvertedRepeat(homologs[0]["intergenic"])[0]
        print(operators)
        #print("inverted repeat found: "+ operator)
        
        #for i in homologs:
        #    i["operator"], i["score"] = findOperatorInIntergenic(i["intergenic"], operator)

        #consensus_data = getConsensus(homologs)
        #print(consensus_data)

    
    #invRepeat = getBestInvertedRepeat(homologs)
    #print(invRepeat)
