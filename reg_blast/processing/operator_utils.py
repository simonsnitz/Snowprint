import pickle
from Bio.pairwise2 import align
#from pprint import pprint
from pathlib import Path

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

def findOperatorInIntergenic(intergenic, operator, ext_length=0):

    operator_length = len(operator)
    
        # This function is needed to extract the ENTIRE aligned region.
        # Otherwise, mismatched ends will not be extracted.
    
    def extractOperator(intergenic, op_align, ext_length):
        begin = 0
        for i in op_align:
            # Find starting position of operator within intergenic region
            if i == '-':
                begin += 1
            else:
                break
        end = (begin + operator_length)
            # Can change indexing of output sequence to include more or less of alignment
        try:
            upstream = intergenic[begin-ext_length:begin].lower()
            mid = intergenic[begin:end]
            downstream = intergenic[end:end+ext_length].lower()
            operator = upstream+mid+downstream
        except:
            operator = intergenic[begin:end]
        return operator


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
        operator = extractOperator(upstr_align, op_align, ext_length)
        return {"operator":operator, "score":score}
    else:
        print('alignment score is not above cutoff')
        return {"operator":None,"score":0}





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
                    allScores.append({"seq":seq,"score":score})

                except:
                    allScores.append({"seq":"None","score":0})
    
    best_score = max([i["score"] for i in allScores])            
    best_operators = [i for i in allScores if i["score"] == best_score]
    
    return best_operators


def findBestPalindrome(intergenic, shortest, longest, winScore, lossScore):

    operators = []
    intergenic = intergenic.upper()
    
    for i in range(shortest, longest):
        ops = findImperfectPalindromes(intergenic, i, winScore, lossScore)
        for j in ops:
            operators.append(j)

    max_score = max([i["score"] for i in operators])
    best_operators = [i for i in operators if i["score"] == max_score]

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

    return {"motif_data":consensus_data, "num_seqs":num_seqs}



def get_consensus_score(operator, consensus_data):
    
    max_score = 0
    consensus_score = 0
    for i in range(0,len(operator)):
        if operator[i].isupper():
            max_score += 1
            consensus_score += consensus_data["motif_data"][i]['score']**2
        else:
            pass
    
    score = round((consensus_score/max_score)*100, 3)

    return score



'''
Input: File with homologs, alignment method

Output: Array with consensus data for each input operator.

'''

def appendOperatorMetadata(acc, to_align, perc_ident):
    
    # Path for operator data
    p = Path('./cache/operator_data/')
    filename = acc+".pkl"
    fp = p / filename

    if fp.is_file():
        print("operator data is already cached for "+str(acc))
        with open(fp, mode="rb") as f:
            operator_data = pickle.load(f)
            return operator_data
    else:

        homologListFile = f"cache/homolog_metadata/updated_metadata/{acc}.pkl"

        with open(homologListFile, mode="rb") as f:
            homologList = pickle.load(f)
            
            if to_align == "find_inverted_repeat":

                    # Iterates this function through multiple relevant scoring parameters
                operators = []

                test_params = [{"w":2,"l":-2}, {"w":2,"l":-3}, {"w":2,"l":-4}]

                for i in test_params:
                    ops = [findBestPalindrome(intergenic=homologList[0]["intergenic"], \
                    shortest=5, longest=15, winScore=i["w"], lossScore=i["l"])][0]
                    for operator in ops:
                        operators.append(operator)
                

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
                "accession": str(acc),
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
                    try:
                        homolog = {}
                        homolog["operator"] =  findOperatorInIntergenic(i["intergenic"], \
                            operator["seq"])["operator"]
                        homolog["score"] =  findOperatorInIntergenic(i["intergenic"], \
                            operator["seq"])["score"]
                        homolog["identity"] = i["identity"]

                        metrics.append(homolog)
                    except:
                        print('error parsing homologList')

                consensus = getConsensus(metrics, alignment_type=to_align)

                consensus_score = get_consensus_score(operator["seq"], consensus)
                #print(operator["seq"], consensus_score, consensus["num_seqs"])
                
                operator["seq"] = findOperatorInIntergenic(homologList[0]["intergenic"], \
                    operator["seq"], 3)["operator"]

                    # Warning: Only the CONSENSUS SCORE is used to identify the 
                        # best operator. This should also incorporate the
                        # NUMBER OF ALIGNMENTS as a metric to make this decision.

                if consensus_score > operator_data["consensus_score"]:
                    operator_data["consensus_score"] = consensus_score
                    operator_data["input_seq"] = operator["seq"]
                    operator_data["num_seqs"] = consensus["num_seqs"]
                    operator_data["motif"] = consensus["motif_data"]
                    
            # Cache operator data
            #p = Path('./cache/operator_data/')
            #filename = acc+".pkl"
            #fp = p / filename

            with fp.open("wb") as f:
                pickle.dump(operator_data, f)
                print("operator data for "+acc+" cached.")
            
            return operator_data


            # Returned as a list for compatability with downstream display script.
            # I programatically choose the best operator, thus, this is not needed.
            # Adapt the downstream script and return this as a dictionary, not a list.
    
    #with open(f'cache/operator_data/{acc}.pkl', mode="wb") as f:
    #    homologList = pickle.load(f)
    #    return [operator_data]       

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
