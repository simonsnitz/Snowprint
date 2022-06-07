import sys
import argparse

fasta = sys.argv[1]



def count_fastas(fasta):
    with open(str(fasta), "r") as f:
        c = 0
        for line in f:
            if line[0] == ">":
                c += 1
        print("There are "+str(c)+" total fasta entries")


if __name__ == "__main__":
    count_fastas(fasta)
