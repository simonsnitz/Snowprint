import pickle
from Bio.pairwise2 import align
from pprint import pprint


"""

Function for finding a homologous operator within a given intergenic region 
via alignment with a known operator

        Do a local alignment with gap penalties:
            Identical characters add 2 points
            Non-Identical characters deduct 0.5 point
            Opening Gap deduct 100 points
            Extending a Gap  deduct 0 points
        Arguments:
            intergenic: Intergenic DNA region (string)
            operator: Known operator for regulator or its homolog (string)

TODO:
    score_cutoff is arbitrarily set to 10% of max score. Is that relevant?
"""

def findOperatorInIntergenic(intergenic, operator):

    operator_length = len(operator)
    
        # This function is needed to extract the ENTIRE aligned region.
        # Otherwise, mismatched ends will not be extracted.
    
    def extractOperator(intergenic, operator):
        beginning = 0
        for i in operator:
            if i == '-':
                beginning += 1
            else:
                break
        end = (beginning + operator_length)
            # Can change position of output sequence to include more or less
        return intergenic[beginning:end]


    try:
        upstr_align, op_align, score, startPos, endPos = \
        align.localms(intergenic, operator, 2, -0.5, -100, 0)[0]
    except:
        return "no intergenic region found"
    
        # Heavily penalizing gap opening avoids errors with downstream data processing, but may miss out on interesting biological features
        # Returns the aligned operator sequence if a similarity threshold is met. Score threshold (7) should be tuned.
        # Set score cutoff to be 10% of max. Arbitrary, but seems reasonable.
    max_score = 2*operator_length
    score_cutoff = max_score*0.1

    if score > score_cutoff:
        operator = extractOperator(upstr_align, op_align)
        return [operator, score]
    else:
        print('alignment score is not above cutoff')
        return [None,0]





    # Functions for finding inverted repeats within a given sequence
'''
The two functions (findImperfectPalindromes and findBestPalindrome) are used to find
    inverted repeats more accurately than the previous version
    (findInvertedRepeat and getBestInvertedRepeat) by analyzing inverted repeats using a 
    scoring function. This allows for finding IRs with mismatches. Also, it can variably 
    penalize IRs with unrealistically long spacer sequences. All reward/penalty scoring
    values can be modified.

TODO:
    - Test limits of intergenic length. A ~300 bp intergenic region finished in ~ 5 seconds.
        I believe I have a hard limit of 500 bp for max intergenic length (getIntergenic.py)
        so this shouldn't be a problem.
    - Optimize the scoring parameters to consistently identify experimentally validated operators.
'''

def complement(sequence):
    compDNA = {"A":"T", "C":"G", "T":"A", "G":"C"}
    complement = ""
    for i in sequence.upper():
        try:
            complement += compDNA[i]
        except:
            print('non standard base found when running complement() function')
            break
    return complement


def findImperfectPalindromes(intergenic, size, winScore, lossScore):

    spacer_penalty = {0:4, 1:4, 2:4, 3:4, 4:4, 5:2, 6:2, 7:0, 8:0, 9:-2, 10:-2, \
        11:-4, 12:-4, 13:-6, 14:-6, 15:-8, 16:-8, 17:-10, 18:-10, 19:-12, 20:-12}

    IRs = (len(intergenic)+1)-(2*size)
    allScores = []
    for i in range(0,IRs):
        repeat = intergenic[i:i+size]
        for j in range(0,IRs):
            # j represents the size of the spacer between inverted repeats
            if j < 21:
                try:
                    compare = complement(intergenic[i+size+j:i+(2*size)+j])[::-1]
                    score = 0
                    for k in range(0,len(repeat)):
                        if repeat[k] == compare[k]:
                            score += winScore
                        else:
                            score += lossScore

                    score += spacer_penalty[j]
                    seq = repeat + intergenic[i+size:i+size+j].lower() + complement(compare)[::-1]
                    allScores.append([seq,score])

                except:
                    allScores.append([0,0,0])
    
    best_score = max([i[1] for i in allScores])            
    best_operators = [i for i in allScores if i[1] == best_score]
    
    return best_operators


def findBestPalindrome(intergenic, shortest, longest, winScore, lossScore):

    operators = []
    intergenic = intergenic.upper()
    
    for i in range(shortest, longest):
        ops = findImperfectPalindromes(intergenic, i, winScore, lossScore)
        for j in ops:
            operators.append(j)

    max_score = max([i[1] for i in operators])
    best_operators = [i for i in operators if i[1] == max_score]

    return best_operators




''' 
Input: homolog operators

Output: dictionary with base/frequency for each position.

'''

def getConsensus(metrics, alignment_type, max_ident=100, min_ident=50):      #DUUUDE! Mad list comprehension!!!
    
        # Filter by identity and by alignment score. Some alignments from ~80% homologs have crap scores
    allOperators = [ i["operator"] for i in metrics
            if i["identity"] >= min_ident and i["identity"] <= max_ident and i["score"] != 0]

    num_seqs = len(allOperators)
	    # Initialize list of dictionaries
    baep = [{base:1}
            for base in allOperators[0] 
        ]

	# Populate dataset with base representation from all input operators
    for operator in allOperators[1:]:
        if len(operator) == len(allOperators[0]):
            for pos in range(0, len(operator)):
                base = operator[pos]

                try:
                    baep[pos][base] +=1
                except:
                    baep[pos].update({base:1})



    def trim_interoperon(baep):
            # Count the number of times a base was not aligned for each position
        spaces = [i['-'] if "-" in i else 0 for i in baep]

        min_num_spaces = min(spaces)
            # Find the start and stop indices for when aligned spaces is minimal
        start = spaces.index(min_num_spaces)
        end = len(spaces) - spaces[::-1].index(min_num_spaces) - 1
        
        return baep[start:end]

    if alignment_type == "whole_interoperon":
        baep = trim_interoperon(baep)



    max_values = [max(baep[pos].values()) for pos in range(0,len(baep))]
    max_score = max(max_values)
	    # Convert base conservation scores as a percent of max
    max_values_percent = [ round(i/max_score,2) for i in max_values ] 

    def get_key(my_dict,val):
        for key, value in my_dict.items():
            if val == value:
                return key

        return "key doesn't exist"

	    # Create a list of most conserved bases at each position
    consensusSeq = [ get_key(baep[pos], max_values[pos])
        for pos in range(0, len(baep))
    ]

	    # Dictionary containing the base and it's score at each position
    consensus_data = [{"base":consensusSeq[i] , "score":max_values_percent[i]} 
            for i in range(0,len(max_values))
        ]

    return [consensus_data, num_seqs]



def get_consensus_score(operator, consensus_data):
    
    max_score = 0
    consensus_score = 0
    for i in range(0,len(operator)):
        if operator[i].isupper():
            max_score += 1
            consensus_score += consensus_data[0][i]['score']**2
        else:
            pass
    
    score = (consensus_score/max_score)*100

    return score



'''
Input: File with homologs, alignment method

Output: Array with consensus data for each input operator.

'''

def appendOperatorMetadata(homologListFile, to_align, perc_ident):
    
    with open(f'{homologListFile}', mode="rb") as f:
        homologList = pickle.load(f)
        
        if to_align == "find_inverted_repeat":

                # TODO:
                # Iterate this function through multiple relevant scoring parameters
            operators = [findBestPalindrome(intergenic=homologList[0]["intergenic"], \
                shortest=5, longest=15, winScore=2, lossScore=-3)][0]

        elif to_align == "whole_interoperon":
            operators = [homologList[0]["intergenic"]]
        else:
            operators = [to_align]
        

        homologList = [i for i in homologList if i["identity"] > perc_ident]

            # Identify the homolog with the lowest percent identity
            # This should instead be a list of all identities. Then I can better 
                # evaluate how BALANCED the identity distribution is (better metric).
        lowest_identity = str(homologList[-1]["identity"])
			
            # Output data to be returned 
        operator_data = { 
            "input_seq": "None",
            "lowest_perc_ident":lowest_identity,
            "num_seqs": "None",
            "consensus_score": 0,
            "motif": "None"
        }
            # Iterate through predicted operators. Pick one with best consensus score
        for operator in operators:
            metrics = []
            for i in homologList:
                
                homolog = {}
                homolog["operator"] =  findOperatorInIntergenic(i["intergenic"], \
                    operator[0])[0]
                homolog["score"] =  findOperatorInIntergenic(i["intergenic"], \
                    operator[0])[1]
                homolog["identity"] = i["identity"]

                metrics.append(homolog)

            consensus = getConsensus(metrics, alignment_type=to_align)

            consensus_score = get_consensus_score(operator[0], consensus)
            print(operator[0], consensus_score, consensus[1])

                # Warning: Only the CONSENSUS SCORE is used to identify the 
                    # best operator. This should also incorporate the
                    # NUMBER OF ALIGNMENTS as a metric to make this decision.

            if consensus_score > operator_data["consensus_score"]:
                operator_data["consensus_score"] = consensus_score
                operator_data["input_seq"] = operator[0]
                operator_data["num_seqs"] = consensus[1]
                operator_data["motif"] = consensus[0]
                  
        pprint(operator_data)

            # Returned as a list for compatability with downstream display script.
            # I programatically choose the best operator, thus, this is not needed.
            # Adapt the downstream script and return this as a dictionary, not a list.
        return [operator_data]       

# What makes a good operator?
    # 1) High consensus score
    # 2) High number of aligned sequences
    # 3) Low / balanced sequence similarity between homologs
        # How do I best evaluate this last metric?

# These should all be displayed to the user
# But there should also be an aggregate quality score
    # A function should be written for each metric 
    # to best normalize them prior to aggregation
        # for example. 40 alignments is likely just as good as 50
        # but 20 alignments is much better than 10


    

if __name__ == "__main__":

    ramr_operator = "ATGAGTGAGTACGCAACTCAT"
    camr_operator = "TATATCgcaGATATA"
    ttgr_operator = "gggggcaATAtTTACAaaCAaccaTGaaTGTAAgTATcacaagatacggggg"
    edcr_intergenic = "GATCCGTTTTCTGCCCGTTGCATGCTATTGTCGCAACGGCTGCTTGAACAAGACAAGGGGATGCGCAACGGGTAGTACCCATGGTCTGGCGCTTCATGGCTGCTAAAGCGCATTTTGGCTGCTATCGTTGGATGGCGGTGAATCAGACAGGCGCGCATGACAGGACCACGTCACAGTTTTCGGATTGCGTGTCGCGATATTTACAGTACAATTGTTACCGTAAAATGCGGCGCGCCCGAAAGTGGGTAAGAACGCAGTATCACGATGGCTTGTTGTTGCCCGATGCCCTGCCTAACCC"

    # print([findBestPalindrome(intergenic=edcr_intergenic, shortest=5, longest=15, \
    #     winScore=2, lossScore=-2)[0][0]])


    bm3r1 = "WP_013083972.pkl"
    test = "WP_011029378.pkl"
    camr = "BAA03510.pkl"
    acc = bm3r1
    
    metadata = appendOperatorMetadata(f'../cache/homolog_metadata/{acc}', \
        'find_inverted_repeat', 50)
