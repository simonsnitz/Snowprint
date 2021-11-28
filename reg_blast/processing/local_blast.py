from Bio.Blast.Applications import NcbiblastpCommandline
import time

# input fasta file
fasta_input = "../cache/blast_DBs/ramr.fasta"

# internal tetr cluster blast database
blast_db = "../cache/blast_DBs/TetR_clusters"


# run blast
start = time.time()
results = NcbiblastpCommandline(query=fasta_input, db=blast_db, outfmt="6 sseqid pident qcovs", num_alignments=1)()[0]
end = time.time()


print("local BLAST took "+str(end-start)+" seconds")

rep_protein = results.split("|")[1]

identity = results.split("|")[2].split("\t")[1]

coverage = results.split("|")[2].split("\t")[2].strip()

print("representative protein: "+rep_protein)

print("identity: "+identity)

print("coverage: "+coverage)