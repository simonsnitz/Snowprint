import sys

fasta = sys.argv[1]

with open(str(fasta), "r") as f:
    c = 0
    for line in f:
        if line[0] == ">":
            c += 1

    print(c)