import requests
import json
from pathlib import Path
import time
import re
from pprint import pprint

from sqlalchemy import create_engine, MetaData, insert
from sqlalchemy.orm import sessionmaker

from src.definitions.define_operon import getOperon


headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.163 Safari/537.36'} 


    # file locations:
tmp = Path("./cache/tmp")
genomes_tmp = tmp / "genomes.txt"




    # returns a 40K fragment of the genome surrounding the target
def acc2genome_frag(genome_id, startPos, stopPos):

        # Define the genome fragment size. Emperically determined based on max operon lengths in the DB
    startPos = startPos - 20000
    if startPos < 0:
        startPos = 0

    stopPos = stopPos + 20000
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore"
    response = requests.get(base_url+"&id="+str(genome_id)+"&seq_start="+str(startPos)+"&seq_stop="+str(stopPos)+"&rettype=fasta_cds_aa")

    if response.ok:
        data = response.text
        with open(genomes_tmp, mode='w+') as f:
            f.write(data)
    else:
        # print(response.status_code)
        print('FATAL: genomes efetch failed. Re-trying ...')
        return

        # Cache genomes text file
    with open(genomes_tmp, mode='r') as f:
        genome_frag = f.read().split("\n")[:-2]

    return genome_frag
    




    # Pull out the gene metadata and the location of the target regulator
        # NOTE: Tried using a 3rd regex with the accession ID, which threw an error.
def parseGenome(genome_frag: str, start: str, stop: str):
    re1 = re.compile(start)
    re2 = re.compile(stop)
    geneIndex = 0
    allGenes = []
    for i in genome_frag:
        if i[0] == '>':
            if re1.search(i):
                if re2.search(i):
                    regIndex = geneIndex
            geneIndex += 1
            allGenes.append(i)
    try:
        return allGenes, regIndex
    except:
        print('FATAL: regulator not found in genome\n')
        return None, None
                



    # Main function that creates Operon and Association records if they don't already exist
def create_operons(acc: str):

        # Point to the sqlite database
    engine = create_engine('sqlite:///cache/Snowprint.db')

        # Connect to db, create session, and access the Alignment table
    conn = engine.connect()
    Session = sessionmaker(bind=engine)
    session = Session()
    meta_data = MetaData(bind=engine)
    MetaData.reflect(meta_data)
    Regulator = meta_data.tables["regulator"]
    Operon = meta_data.tables["operon"]
    Association = meta_data.tables["association"]
    Alignment = meta_data.tables["alignment"]


            # Extract Alignment record associated with accession ID argument
    record = session.query(Alignment).filter_by(query_id=acc).first()
    if record == None:    
        print('FATAL: No alignment found for '+str(acc))
        return

        # Extract accession IDs for homologs within the alignment
    else:
        homologs = json.loads(record.homologs)
        accessions = [i['accession'] for i in homologs]

        # Pull out regulators associated with the alignment
    regulators = [session.query(Regulator).filter_by(prot_id=acc).first() for acc in accessions]
        # Filter out empty Regulator records
    regulators = [reg for reg in regulators if reg != None]


        # Get operon data for each regulator if it doesn't already exist
    for reg in regulators:
        assoc = session.query(Association).filter_by(regulator_id=reg.id).first()
        if assoc == None:

            genome_frag = acc2genome_frag(reg.genome_id, reg.start_pos, reg.stop_pos)

            allGenes, index = parseGenome(genome_frag, str(reg.start_pos), str(reg.stop_pos))

            if allGenes != None:

                operon_genes, regIndex = getOperon(allGenes, index, reg.start_pos, reg.strand)
                
                operon_start = operon_genes[0]["start"]
                operon_stop = operon_genes[-1]["stop"]

                operon = {"prot_id": reg.prot_id, "genome_id": reg.genome_id, \
                     "operon": operon_genes, "regIndex": regIndex, "start": operon_start,\
                        "stop": operon_stop}

            else:
                print('FATAL: Genome parsing failed')
                return

                # Check to see if a DB record already exists for the operon.
            operon_record = session.query(Operon).filter_by(
                genome_id= operon["genome_id"],
                start_pos= operon["start"],
                stop_pos= operon["stop"]).first()

            if operon_record == None:

                    # Create a new operon record
                new_operon = (
                    insert(Operon).values(
                    operon= json.dumps(operon["operon"]),
                    start_pos= operon["start"],
                    stop_pos= operon["stop"],
                    genome_id= operon["genome_id"]
                        )
                )
                conn.execute(new_operon)
                # print('UPDATE: Added an operon entry for '+str(operon["prot_id"]))

            else:
                continue

                # Pull out regulator and operon records to associate with each other
            reg_to_associate = session.query(Regulator).filter_by(prot_id=operon["prot_id"]).first()
            operon_to_associate = session.query(Operon).filter_by(operon=json.dumps(operon["operon"])).first()

                # Check to see if a DB record already exists for the association.
                    # There may be a case where an operon exists, but a reg-operon association doesn't
            assoc_record = session.query(Association).filter_by(
                regulator_id= reg_to_associate.id,
                operon_id= operon_to_associate.id).first()


            if assoc_record == None:

                    # Create an new association record
                new_association = (
                    insert(Association).values(
                    regulator_id = reg_to_associate.id,
                    operon_id = operon_to_associate.id,
                    reg_index= operon["regIndex"]
                        )
                )
                new_association.regulator = reg_to_associate
                new_association.operon = operon_to_associate

                conn.execute(new_association)

                # print('UPDATE: Added an association entry for '+str(operon["prot_id"]))
            
            else:
                continue
    
    conn.close()

        #TODO:
            # find out how to add to an association table, properly ...



if __name__ == "__main__":

    #fetch_operon_data("WP_000113282.1")

            # Point to the sqlite database
    engine = create_engine('sqlite:///cache/Snowprint.db')

        # Connect to db, create session, and access the Alignment table
    conn = engine.connect()
    Session = sessionmaker(bind=engine)
    s = Session()
    meta_data = MetaData(bind=engine)
    MetaData.reflect(meta_data)

    Operon = meta_data.tables["operon"]

    operons = s.query(Operon)

    op_length = [op.stop_pos - op.start_pos for op in operons]
    print(len(op_length))
    op_filt = [op for op in op_length if op < 16000]
    print(len(op_filt))
    print(len(op_filt)/len(op_length))