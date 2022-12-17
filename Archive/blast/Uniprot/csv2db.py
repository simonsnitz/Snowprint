import csv
from re import I

def create_fasta_db():
    with open('uniprot-tetr+regulator.tab', 'r') as f:
        reader = csv.reader(f)
        c = 0
        with open('uniprot_TetRs.fasta', "w+") as out:
            entries = ""
            for row in reader:
                if c != 0:
                    data = row[0].split('\t')
                    fasta = ">"+str(data[0])+"\n"
                    sequence = str(data[5])+"\n"
                    entries += fasta
                    entries += sequence
                c += 1
            out.write(entries)


def create_fasta_map():
    with open('uniprot-tetr+regulator.tab', 'r') as f:
        reader = csv.reader(f)
        c = 0
        with open('uniprot_TetRs_taxid_map.txt', "w+") as out:
            entries = ""
            for row in reader:
                if c != 0:
                    data = row[0].split('\t')
                    fasta = str(data[0])+"\t"
                    sequence = str(data[3])+"\n"
                    entries += fasta
                    entries += sequence
                c += 1
            out.write(entries)



def check_fasta():
    with open("uniprot_TetRs.fasta", "r") as f:
        c = 0
        for line in f:
            if line[0] == ">":
                c += 1

        print(c)


if __name__ == "__main__":
    create_fasta_map()