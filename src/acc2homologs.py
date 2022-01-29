from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.NCBIXML import read, parse

from sqlalchemy import create_engine, MetaData, insert
from sqlalchemy.orm import sessionmaker

import requests
import json
from pathlib import Path

    # Cache locations:
cache = Path("./cache")
blast_tmp = cache / "blast_results.xml"
alignment_tmp = cache / "alignment.json"



def accID2sequence(accID):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("bad request "+ str(response.status_code))



    # Hitlist_size trumps perc_ident
def blast_and_cache(sequence, hitlist_size=100, db="nr"):
        with open(f'cache/blast_results.xml', mode="w+") as f:
        #with open(blast_tmp, mode="w+") as f:      # did not work. Why???
            print("Starting BLAST")
            blast_results = qblast("blastp", db, sequence, hitlist_size=hitlist_size)
            f.write(blast_results.read())
            print('cached blast result')
            return blast_results.read()



def homologs2metadata():

    with open(blast_tmp, mode="r") as f:

        records = [record for record in parse(f)]
        for record in records:
            homologs = [
                {"accession":alignment.accession, 
                "identity":round((hsp.identities/hsp.align_length)*100, 2),
                "coverage": round((hsp.align_length/record.query_length)*100, 2)
                }
                for alignment in record.alignments
                    for hsp in alignment.hsps
            ]

    alignment = json.dumps(homologs, indent=4)
    print(alignment)
  
    with open(alignment_tmp, 'w+') as f:
        f.write(alignment)
        print('homolog acc/pident/qcov cached')
    
    return alignment



def add_alignment_to_db(acc):
    
        # Create alignment data from blast cache
    alignment = homologs2metadata()

        # Point to sqlite database
    engine = create_engine('sqlite:///API/GroovIO.db')

        # Connect to db and pull Alignment table
    conn = engine.connect()
    meta_data = MetaData(bind=engine)
    MetaData.reflect(meta_data)
    Alignment = meta_data.tables["alignment"]

        # Create a new record with accession ID argument and alignment data
    new_row = (
        insert(Alignment).values(
            query_id=acc,
            homologs = alignment)
        )
    
        # Add the new record and commit it to the DB
    result = conn.execute(new_row)
    conn.close()
    


def query_db_alignment(acc):

        # Point to sqlite database
    engine = create_engine('sqlite:///API/GroovIO.db')

        # Connect to db. create session and pull Alignment table
    conn = engine.connect()
    Session = sessionmaker(bind=engine)
    session = Session()
    meta_data = MetaData(bind=engine)
    MetaData.reflect(meta_data)
    Alignment = meta_data.tables["alignment"]

        # Extract record associated with accession ID argument
    record = session.query(Alignment).filter_by(query_id=acc).first()

        # Do something with extracted record
    print(json.loads(record.homologs)[40])

    conn.close()



if __name__=="__main__":

        #alkx
    acc = "AEM66515"


