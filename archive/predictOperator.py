



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

#def getOperator(intergenic, size):
def intergenic2Operators(intergenic, size):
    
    operators = []

    for i in range(0,len(intergenic)-((2*size)+5)):
        seq = intergenic[i:i+size]
        revCompSeq = complement(seq)[::-1]
        maxAllowedSpacer = i+size+5
            
        for j in range(i+size,maxAllowedSpacer):
            compare = intergenic[j:j+size]
            if compare == revCompSeq:
                
                seqF = seq.upper()
                middle = intergenic[i+size:j].lower()
                seqR = compare.upper()
                operator = seqF+middle+seqR
                operators.append(operator)
    return operators

def createOperatorList(intergenic):
    stringency = [9,8,7,6,5,4]
    if intergenic != None:
        for i in stringency:
            operators = intergenic2Operators(intergenic, i)
            if operators:
                return operators
    else:
        return ""


if __name__ == "__main__":

    with open("intergenic.txt") as f:
        intergenic = ""
        for i in f.readlines():
            intergenic += i

    print(intergenic)
    operators = createOperatorList(intergenic)
    print(operators)
