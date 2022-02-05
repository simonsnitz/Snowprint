from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.NCBIXML import parse

from sqlalchemy import create_engine, MetaData, insert
from sqlalchemy.orm import sessionmaker

from Bio.Blast.Applications import NcbiblastpCommandline

import requests
import json
import time
from pathlib import Path



    # Cache locations:
cache = Path("./cache")
blast_tmp = cache / "blast_results.xml"
alignment_tmp = cache / "alignment.json"


    # Input protein accession ID, output sequence in fasta format
def accID2sequence(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("bad request "+ str(response.status_code))


    # Input protein sequence. Output cached blast results
def blast(acc: str, num_aligns=100):

    print("Starting BLAST")


    blast_db = "blast/TetR/TetRs"

    seq = accID2sequence(acc)

    blast_cline = NcbiblastpCommandline(db=blast_db, outfmt="6 sseqid pident qcovs", num_alignments=num_aligns)


    results, err = blast_cline(stdin=seq)

    results = results.split("\n")[:-1]
    
    homologs = [{"accession": r.split("|")[1], \
            "identity": r.split("|")[2].split("\t")[1], \
            "coverage": r.split("|")[2].split("\t")[2].strip()} \
             for r in results ]

    return homologs

    # with open(f'cache/blast_results.xml', mode="w+") as f:
    #with open(blast_tmp, mode="w+") as f:      # did not work. Why???

            # blast_results = qblast("blastp", db, sequence, hitlist_size=hitlist_size)
            # f.write(blast_results.read())
            # print('cached blast result')
            # return blast_results.read()


    # Input blast alignment results, output JSON dict with homolog data
def homologs2metadata():

    with open(blast_tmp, mode="r") as f:
        try:
            records = [record for record in parse(f)]
            for record in records:
                homologs = [
                    {
                        "accession":alignment.accession, 
                        "identity":round((hsp.identities/hsp.align_length)*100, 2),
                        "coverage": round((hsp.align_length/record.query_length)*100, 2)
                    }
                    for alignment in record.alignments
                        for hsp in alignment.hsps
                ]

            alignment = json.dumps(homologs, indent=4)
            return alignment

        except IndexError:
            print('There is an issue with the cached BLAST record')



    # Checks if Alignment record already exists in DB.
        # If not, fetches the sequence, BLASTs it, creates an alignment, and adds it to the DB
def create_alignment(acc: str):

        # Point to sqlite database
    engine = create_engine('sqlite:///API/GroovIO.db')

        # Connect to db, create session, and access the Alignment table
    conn = engine.connect()
    Session = sessionmaker(bind=engine)
    session = Session()
    meta_data = MetaData(bind=engine)
    MetaData.reflect(meta_data)
    Alignment = meta_data.tables["alignment"]


        # Extract record associated with accession ID argument
    record = session.query(Alignment).filter_by(query_id=acc).first()

    if record == None:    
        print('Alignment not found for '+str(acc))
            # Fetch sequence, blast it, and create alignment data
        alignment = blast(acc)
        #alignment = homologs2metadata()
            # Create a new record with accession ID and alignment data
        new_row = (
            insert(Alignment).values(
                query_id=acc,
                homologs = alignment
                )
        )
            # Add the new record and commit it to the DB
        # conn.execute(new_row)
        print('Added an alignement for '+str(acc)+' to the DB')

    else:
        print('DB alignment already exists for '+str(acc))

    conn.close()

#TODO:
    # find a way to check if blast cache is for the query regulator



if __name__=="__main__":

    blast_db = "../blast/TetR/TetRs"

    seq = accID2sequence("WP_000113282.1")

    start = time.time()
    blast_cline = NcbiblastpCommandline(db=blast_db, outfmt="6 sseqid pident qcovs", num_alignments=100)


    results, err = blast_cline(stdin=seq)

    results = results.split("\n")[:-1]
    homologs = [{"protein_id": r.split("|")[1], \
            "identity": r.split("|")[2].split("\t")[1], \
            "coverage": r.split("|")[2].split("\t")[2].strip()} \
             for r in results ]

    #rep_protein = results.split("|")[1]
    #identity = results.split("|")[2].split("\t")[1]
    #coverage = results.split("|")[2].split("\t")[2].strip()

    stop = time.time()

    print(results)
    print(homologs[0]["protein_id"])
    #print("took "+str(stop-start)+" seconds")