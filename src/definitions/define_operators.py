

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
    
    if len(allScores) > 0:      # sometimes you get an empty array
        best_score = max([i["score"] for i in allScores])            
        best_operators = [i for i in allScores if i["score"] == best_score]
        
        return best_operators





def findBestPalindrome(intergenic, shortest, longest, winScore, lossScore):

    operators = []
    intergenic = intergenic.upper()
    
    for i in range(shortest, longest):
        ops = findImperfectPalindromes(intergenic, i, winScore, lossScore)
        if ops != None:
            for j in ops:
                operators.append(j)

    max_score = max([i["score"] for i in operators])
    best_operators = [i for i in operators if i["score"] == max_score]

    return best_operators